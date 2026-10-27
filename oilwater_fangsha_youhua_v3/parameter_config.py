#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parameter_config.py —— Python 脚本读取 config/*.env 参数的公共模块
======================================================================

【为什么需要这个文件】
  原交付目录里 analyze_production.py、audit_parallel_bond.py、make_water_fangsha_pack.py、
  validate_water_fangsha_pack.py 都写了 `from parameter_config import ...`，但这个文件缺失，
  所以 run_production.sh 在"几何校验"一步就会报 ModuleNotFoundError（审计 P1-6）。
  这里按 C++ 端 parameterConfig.h 的规则把它补上，保证三方（bash / C++ / Python）读到同一套值。

【读取规则】（与 parameterConfig.h、config/load_parameters.sh 一致）
  1. 先读 PARAMETER_DEFAULTS_FILE（默认 config/defaults.env）——它是全部合法参数名的清单；
  2. 再读 CONFIG_FILE（实验剖面）。剖面里可以用 BASE_CONFIG=相对路径 继承另一个剖面，
     继承链从最底层基类开始依次覆盖；
  3. 同名环境变量优先级最高（例如命令行 `V_INF=0.06 bash run_production.sh`）。
  4. 行格式 KEY=VALUE，# 开头为注释，值两端的引号会被去掉；未知参数名直接报错（防拼写错误）。

【用法】
  from parameter_config import load_parameters, real, text, integer
  cfg = load_parameters()
  lx = real(cfg, "DOMAIN_LX_M")      # 转成 float
  n  = integer(cfg, "N")             # 转成 int
  pf = text(cfg, "PACK_FILE")        # 原始字符串
"""
from __future__ import annotations

import os
import re
from pathlib import Path

_LINE = re.compile(r"^([A-Z][A-Z0-9_]*)=(.*)$")


def _trim(value: str) -> str:
    """去掉两端空白和成对引号。"""
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        value = value[1:-1]
    return value


def _read_file(path: Path) -> dict[str, str]:
    """读取一个 KEY=VALUE 文件，返回字典（不做继承）。"""
    if not path.is_file():
        raise FileNotFoundError(f"parameter file not found: {path}")
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        m = _LINE.match(line)
        if not m:
            raise ValueError(f"invalid parameter line in {path}: {raw}")
        values[m.group(1)] = _trim(m.group(2))
    return values


def _load_recursive(path: Path, allowed: set[str], seen: set[Path], out: dict[str, str]) -> None:
    """按 BASE_CONFIG 继承链加载：先加载基类，再用本文件覆盖。"""
    canonical = path.resolve()
    if canonical in seen:
        raise ValueError(f"cyclic BASE_CONFIG inheritance at {canonical}")
    seen.add(canonical)
    values = _read_file(canonical)
    base = values.get("BASE_CONFIG", "")
    if base:
        base_path = Path(base) if os.path.isabs(base) else canonical.parent / base
        _load_recursive(base_path, allowed, seen, out)
    for key, value in values.items():
        if key not in allowed:
            raise KeyError(f"unknown parameter '{key}' in {canonical}; check spelling against defaults.env")
        out[key] = value


def load_parameters() -> dict[str, str]:
    """返回合并后的参数字典（文件值 + 环境变量覆盖）。"""
    defaults = Path(os.environ.get("PARAMETER_DEFAULTS_FILE") or "config/defaults.env")
    profile = os.environ.get("CONFIG_FILE") or ""
    config = _read_file(defaults)                 # defaults.env 本身定义了合法参数名
    allowed = set(config.keys())
    if profile and Path(profile).resolve() != defaults.resolve():
        _load_recursive(Path(profile), allowed, set(), config)
    for key in list(config.keys()):               # 环境变量覆盖
        env = os.environ.get(key)
        if env:
            config[key] = env
    return config


def text(config: dict[str, str], name: str) -> str:
    if name not in config:
        raise KeyError(f"missing parameter {name}")
    return config[name]


def real(config: dict[str, str], name: str) -> float:
    raw = text(config, name)
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"invalid real parameter {name}={raw}") from exc
    if value != value or value in (float("inf"), float("-inf")):
        raise ValueError(f"invalid real parameter {name}={raw}")
    return value


def integer(config: dict[str, str], name: str) -> int:
    raw = text(config, name)
    try:
        return int(raw, 10)
    except ValueError as exc:
        raise ValueError(f"invalid integer parameter {name}={raw}") from exc


def real_list(config: dict[str, str], name: str) -> list[float]:
    """逗号分隔的实数列表，例如 PV_SAMPLES=1,5,10。"""
    raw = text(config, name)
    items = [float(v) for v in raw.split(",") if v.strip()]
    if not items:
        raise ValueError(f"empty parameter {name}")
    return items


if __name__ == "__main__":          # 直接运行：打印合并后的全部参数，方便检查
    for k, v in sorted(load_parameters().items()):
        print(f"{k}={v}")

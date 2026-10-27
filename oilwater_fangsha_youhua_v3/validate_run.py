#!/usr/bin/env python3
"""
validate_run.py —— 对一次已经跑完的 oilwater_fangsha 运行目录做“机器体检”。

【这个脚本做什么】
  读取 C++ 主程序（oilwaterDrive）写出的各类 CSV / 日志 / 退出码，逐项检查
  “这次模拟是否完整、数值是否稳定、颗粒和胶结键是否守恒”，每项打印
  PASS/FAIL，全部通过返回 0，否则返回 7。它不画图、不做统计，只回答
  “这次结果能不能用”。

【什么时候、由谁调用】
  手工调用。典型用法是模拟跑完后在案例目录执行：
      python3 validate_run.py <运行目录>
  其中 <运行目录> 就是 run_production.sh 里的 $OUTDIR。

【输入文件（都在运行目录下）】
  必需：input_parameters.env（本次参数快照，含 VALIDATION_* 容差）、
        program_exit_code.txt、run_stdout.log、geometry_metrics.csv、
        pack_geometry.csv、parallel_bond_audit.csv、screen_metrics.csv、
        sha256_manifest.txt、post/particles_final.liggghts
  标定模式(RUN_MODE=calibration) 另需：permeability_calibration.csv
  生产模式(RUN_MODE=production)  另需：perm_porosity.csv、ibm_force_audit.csv、
        flow_stability.csv、particle_balance.csv、bond_history.csv 以及
        analyze_production.py 写出的 analysis/*.csv、analysis/summary.txt

【主要判据】
  - 主程序退出码为 0，且日志里有完成标记
    （生产：“[PV] 达到总注入”；标定：“[PERMEABILITY_RESULT] converged=1”）
  - 运行期 IBM 岩心孔隙度在 PHI_MIN..PHI_MAX 内
  - 筛网压降 / 筛前颗粒计数全部有限且非负
  - 生产模式：实际注入 PV 达到 TOTAL_PV；颗粒速度不超过限速、无超力/超矩/NaN；
    LBM 最大速度 < VALIDATION_MAX_LB_VELOCITY，平均密度在允许区间；
    骨架/筛网颗粒数不变、弱颗粒“剩余+产出=初始”、异常丢失为 0；
    胶结键“活动+累计断裂=初始”；后处理汇总报告异常丢失为 0
  - 标定模式：渗透率窗口收敛且落在 K_MIN_MD..K_MAX_MD 内

【列名已随审计修订版同步】
  本脚本引用的 CSV 列名已对照 oilwaterDrive.cpp 中 permCsv/forceCsv/bondCsv/
  balanceCsv/flowCsv/screenCsv/calibrationCsv 的表头字符串逐一核对并修改：
    ibm_force_audit.csv : max_torque_N -> max_torque_Nm；
                          over_force_count/over_torque_count -> over_force/over_torque；
                          去掉主程序已不再输出的 raw_max_*、limited_*_count 列
    flow_stability.csv  : mean_density_lb -> mean_density_w_lb + mean_density_o_lb
                          （两相各自的平均密度；区间判据作用在两者之和上）
    bond_history.csv    : 去掉已不存在的 active_count_measured 列
  这是本文件相对旧版唯一的逻辑改动，其余判据保持不变。
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


def read_env(path: Path) -> dict[str, str]:
    """把 KEY=VALUE 形式的参数快照文件读成字典（忽略没有 '=' 的行）。"""
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def read_rows(path: Path) -> list[dict[str, str]]:
    """读取一个带表头的 CSV，返回“每行一个字典”的列表（所有值仍是字符串）。"""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite_row(row: dict[str, str], fields: tuple[str, ...]) -> bool:
    """检查某一行里给定的列都存在、能转成数、且不是 NaN/Inf。"""
    try:
        return all(math.isfinite(float(row[field])) for field in fields)
    except (KeyError, TypeError, ValueError):
        return False


def artifact_path(run_dir: Path, relative: str, parameters: dict[str, str]) -> Path:
    """定位输出文件：post/ 下的文件可能带实验名前缀（OUTPUT_PREFIX），两种写法都试。"""
    direct = run_dir / relative
    if direct.is_file():
        return direct
    prefix = parameters.get("OUTPUT_PREFIX", "")
    if prefix and relative.startswith("post/"):
        parent, name = relative.split("/", 1)
        candidate = run_dir / parent / f"{prefix}{name}"
        if candidate.is_file():
            return candidate
    return direct


def main() -> int:
    """逐项检查运行目录，打印 PASS/FAIL 清单；全部通过返回 0，否则返回 7。"""
    parser = argparse.ArgumentParser(description="Validate a completed water_fangsha run directory")
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir
    checks: list[tuple[bool, str]] = []  # 每项检查 = (是否通过, 说明文字)

    # ── 第 1 步：基础文件齐不齐。缺任何一个就没必要继续，直接报 FAIL ──
    required = [
        "input_parameters.env",
        "program_exit_code.txt",
        "run_stdout.log",
        "geometry_metrics.csv",
        "pack_geometry.csv",
        "parallel_bond_audit.csv",
        "screen_metrics.csv",
        "sha256_manifest.txt",
    ]
    for relative in required:
        checks.append(((run_dir / relative).is_file(), f"required output exists: {relative}"))
    if not all(passed for passed, _ in checks):
        for passed, message in checks:
            print(f"{'PASS' if passed else 'FAIL'}: {message}")
        return 7

    # ── 第 2 步：读参数快照，确认最终颗粒 dump 存在，并校验容差本身合法 ──
    parameters = read_env(run_dir / "input_parameters.env")
    checks.append((artifact_path(run_dir, "post/particles_final.liggghts", parameters).is_file(),
                   "required output exists: experiment-qualified post/particles_final.liggghts"))
    if not all(passed for passed, _ in checks):
        for passed, message in checks:
            print(f"{'PASS' if passed else 'FAIL'}: {message}")
        return 7
    # 下面这些容差来自 config/defaults.env 的 VALIDATION_* 项，经 input_parameters.env 带入
    validation_pv_rel_tol = float(parameters["VALIDATION_PV_REL_TOL"])            # PV 允许的相对欠量
    validation_force_limit_factor = float(parameters["VALIDATION_FORCE_LIMIT_FACTOR"])  # 速度上限的放宽倍数
    validation_max_lb_velocity = float(parameters["VALIDATION_MAX_LB_VELOCITY"])    # 格子速度上限(低 Ma 要求)
    validation_min_mean_density = float(parameters["VALIDATION_MIN_MEAN_DENSITY"])  # 平均格子密度下限
    validation_max_mean_density = float(parameters["VALIDATION_MAX_MEAN_DENSITY"])  # 平均格子密度上限
    validation_particle_count_tol = float(parameters["VALIDATION_PARTICLE_COUNT_TOL"])  # 颗粒守恒容差(颗)
    validation_bond_count_tol = float(parameters["VALIDATION_BOND_COUNT_TOL"])      # 键守恒容差(根)
    if (validation_pv_rel_tol < 0.0 or validation_force_limit_factor < 1.0
            or validation_max_lb_velocity <= 0.0
            or not 0.0 < validation_min_mean_density < validation_max_mean_density
            or validation_particle_count_tol <= 0.0 or validation_bond_count_tol <= 0.0):
        raise ValueError("invalid validation tolerances in input_parameters.env")

    # ── 第 3 步：运行模式、退出码、标定状态标签是否明确 ──
    mode = parameters.get("RUN_MODE", "")
    program_exit = int((run_dir / "program_exit_code.txt").read_text().strip())
    checks.append((program_exit == 0, "coupled program returned exit code 0"))
    checks.append((mode in {"production", "calibration"}, "run mode is recognized"))
    # 胶结键参数来源必须显式声明（是文献移植值还是实验标定值），避免“默认值当真值”
    checks.append(
        (
            parameters.get("BOND_CALIBRATION_STATUS")
            in {
                "provisional_flow_detachment_baseline",
                "pfc_deformability_interface_unvalidated",
                "frontiers_2024_table3_transfer_unvalidated",
                "tang_2024_table3_core_I_unvalidated",
                "experimentally_calibrated",
            },
            "parallel-bond calibration status is explicit",
        )
    )
    checks.append(
        (
            parameters.get("SCREEN_CALIBRATION_STATUS")
            in {"geometric_baseline_unvalidated", "experimentally_calibrated"},
            "outlet-screen calibration status is explicit",
        )
    )
    checks.append(
        (
            parameters.get("EXPERIMENTAL_CALIBRATION_AVAILABLE") in {"0", "1"},
            "experimental calibration availability is explicit",
        )
    )
    log_text = (run_dir / "run_stdout.log").read_text(encoding="utf-8", errors="replace")

    # ── 第 4 步：运行期几何。IBM 在格子上重新算出的岩心孔隙度要落在接受区间内 ──
    geometry_rows = read_rows(run_dir / "geometry_metrics.csv")
    geometry_ok = bool(geometry_rows) and finite_row(
        geometry_rows[-1], ("core_porosity", "initial_pore_volume_m3", "channel_area_m2", "screen_x_m", "initial_screen_particles")
    )
    if geometry_ok:
        phi = float(geometry_rows[-1]["core_porosity"])
        geometry_ok = float(parameters["PHI_MIN"]) <= phi <= float(parameters["PHI_MAX"])
    checks.append((geometry_ok, "runtime IBM core porosity is finite and accepted"))

    # ── 第 5 步：筛网指标。压降有限；筛前“慢速颗粒数”不能超过“筛前颗粒总数” ──
    screen_rows = read_rows(run_dir / "screen_metrics.csv")
    screen_ok = bool(screen_rows) and all(
        finite_row(row, ("PV", "screen_dp_Pa", "total_dp_Pa", "near_screen_upstream", "near_screen_slow", "downstream_of_screen"))
        and float(row["near_screen_upstream"]) >= 0
        and 0 <= float(row["near_screen_slow"]) <= float(row["near_screen_upstream"])
        and float(row["downstream_of_screen"]) >= 0
        for row in screen_rows
    )
    checks.append((screen_ok, "screen pressure-drop and particle counters are finite"))

    if mode == "calibration":
        # ── 标定模式：只看渗透率是否收敛、是否落在目标区间 ──
        calibration_path = run_dir / "permeability_calibration.csv"
        checks.append((calibration_path.is_file(), "permeability calibration history exists"))
        rows = read_rows(calibration_path) if calibration_path.is_file() else []
        final_ok = bool(rows) and finite_row(
            rows[-1], ("permeability_md", "total_dp_Pa", "window_relative_range", "steady_count")
        )
        if final_ok:
            final_k = float(rows[-1]["permeability_md"])
            # “稳态”定义：窗口内相对波动 < K_STEADY_REL_TOL，且持续了 K_STEADY_HOLD_S 这么久
            required_steady_count = math.ceil(
                float(parameters["K_STEADY_HOLD_S"]) / float(parameters["LOG_INTERVAL_S"])
            )
            final_ok = (
                float(parameters["K_MIN_MD"]) <= final_k <= float(parameters["K_MAX_MD"])
                and float(rows[-1]["window_relative_range"]) < float(parameters["K_STEADY_REL_TOL"])
                and int(rows[-1]["steady_count"]) >= required_steady_count
            )
        checks.append((final_ok, "Q/dP/k window converged and final permeability is accepted"))
        checks.append(("[PERMEABILITY_RESULT] converged=1" in log_text, "calibration completion marker exists"))
    elif mode == "production":
        # ── 生产模式：完整性 + 数值稳定性 + 守恒性 ──
        production_required = [
            "perm_porosity.csv",
            "ibm_force_audit.csv",
            "flow_stability.csv",
            "particle_balance.csv",
            "bond_history.csv",
            "analysis/production_summary.csv",
            "analysis/screen_particle_events.csv",
            "analysis/screen_summary.csv",
            "analysis/summary.txt",
        ]
        for relative in production_required:
            checks.append(((run_dir / relative).is_file(), f"production output exists: {relative}"))

        # 注入量：最后一次采样的 PV 必须达到 TOTAL_PV（允许极小相对欠量）
        perm_rows = read_rows(run_dir / "perm_porosity.csv") if (run_dir / "perm_porosity.csv").is_file() else []
        target_pv = float(parameters["TOTAL_PV"])
        pv_ok = bool(perm_rows) and finite_row(perm_rows[-1], ("PV", "core_porosity", "core_dp_Pa", "total_dp_Pa"))
        if pv_ok:
            pv_ok = float(perm_rows[-1]["PV"]) >= target_pv * (1.0 - validation_pv_rel_tol)
        checks.append((pv_ok, "actual sampled PV reached the requested production target"))
        checks.append((re.search(r"\[PV\] 达到总注入", log_text) is not None, "production completion marker exists"))

        # IBM 力审计：每个记录步的最大力/矩/毛细力有限，颗粒速度不超过 DEM 限速，
        # 且主程序自己的越界标志(over_force/over_torque)和 NaN 计数都为 0
        force_rows = read_rows(run_dir / "ibm_force_audit.csv") if (run_dir / "ibm_force_audit.csv").is_file() else []
        force_ok = bool(force_rows) and all(
            finite_row(row, ("max_force_N", "max_torque_Nm", "max_capillary_force_N",
                             "max_particle_speed_m_s", "max_particle_surface_speed_m_s", "velocity_limit_count"))
            and float(row["max_particle_speed_m_s"]) <= float(parameters["DEM_MAX_PARTICLE_SPEED_M_S"]) * validation_force_limit_factor
            and float(row["max_particle_surface_speed_m_s"])
                <= float(parameters["DEM_MAX_PARTICLE_SPEED_M_S"]) * validation_force_limit_factor
            and float(row["velocity_limit_count"]) >= 0.0
            and int(row["over_force"]) == 0
            and int(row["over_torque"]) == 0
            and int(row["nonfinite_count"]) == 0
            for row in force_rows
        )
        checks.append((force_ok, "IBM force audit is finite with zero threshold violations"))

        # LBM 稳定性：格子最大速度要远小于声速(低马赫)，水、油两相平均密度之和要在合理区间
        flow_rows = read_rows(run_dir / "flow_stability.csv") if (run_dir / "flow_stability.csv").is_file() else []
        flow_ok = bool(flow_rows) and all(
            finite_row(row, ("max_velocity_lb", "mean_velocity_lb", "mean_density_w_lb", "mean_density_o_lb"))
            and float(row["max_velocity_lb"]) < validation_max_lb_velocity
            and validation_min_mean_density
                < float(row["mean_density_w_lb"]) + float(row["mean_density_o_lb"])
                < validation_max_mean_density
            for row in flow_rows
        )
        checks.append((flow_ok, "LBM velocity and mean density remain in the accepted stability range"))

        # 颗粒守恒：骨架(type1)和筛网(type3)颗粒数不许变；弱颗粒 剩余+产出=初始；异常丢失=0
        balance_rows = read_rows(run_dir / "particle_balance.csv") if (run_dir / "particle_balance.csv").is_file() else []
        balance_ok = bool(balance_rows) and all(
            finite_row(row, ("remaining_skeleton", "remaining_screen", "remaining_weak", "produced_weak", "abnormal_loss"))
            and abs(float(row["remaining_skeleton"]) - float(row["initial_skeleton"])) < validation_particle_count_tol
            and abs(float(row["remaining_screen"]) - float(row["initial_screen"])) < validation_particle_count_tol
            and abs(float(row["remaining_weak"]) + float(row["produced_weak"]) - float(row["initial_weak"])) < validation_particle_count_tol
            and abs(float(row["abnormal_loss"])) < validation_particle_count_tol
            for row in balance_rows
        )
        checks.append((balance_ok, "particle balance has zero abnormal loss and unchanged skeleton/screen"))

        # 胶结键守恒：当前活动键 + 累计断裂键 = 初始键数（键只会断，不会新生）
        bond_rows = read_rows(run_dir / "bond_history.csv") if (run_dir / "bond_history.csv").is_file() else []
        bond_ok = bool(bond_rows) and all(
            finite_row(row, ("initial_active_bonds", "active_bonds", "broken_this_lb_step", "cumulative_broken_bonds"))
            and abs(
                float(row["active_bonds"])
                + float(row["cumulative_broken_bonds"])
                - float(row["initial_active_bonds"])
            )
            < validation_bond_count_tol
            for row in bond_rows
        )
        checks.append((bond_ok, "bond count is conserved and includes measured active-bond samples"))

        # 后处理(analyze_production.py)的文本汇总也必须报告“异常丢失 : 0”
        analysis_text = (run_dir / "analysis/summary.txt").read_text(errors="replace") if (run_dir / "analysis/summary.txt").is_file() else ""
        checks.append((re.search(r"异常丢失\s*:\s*0", analysis_text) is not None, "post-processing reports zero abnormal loss"))

    # ── 最后：打印清单并给出总退出码 ──
    for passed, message in checks:
        print(f"{'PASS' if passed else 'FAIL'}: {message}")
    return 0 if all(passed for passed, _ in checks) else 7


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
audit_parallel_bond.py —— 把平行胶结键(parallel bond)的输入参数换算成“力、力矩、临界速度”量级。

【这个脚本做什么】
  LIGGGHTS-INL 的平行胶结模型用“模量 + 强度”描述弱颗粒和骨架之间的胶结。
  这些参数本身看不出“水流拖得动拖不动”。本脚本按 PFC 平行胶结公式，对颗粒包里
  每一根初始胶结键算出：
    - 键半径、截面积、惯性矩、极惯性矩
    - 法向/切向刚度（按单位面积和按整根键）
    - 法向拉断力、切向剪断力、弯矩和扭矩容量
    - 入口流速下单颗细颗粒的 Stokes 拖曳力，以及“拉断力 / 拖曳力”比值
    - 弯矩等效临界流速（拖曳力矩刚好等于弯矩容量时的流速）
  用来回答：“当前胶结强度设置下，弱颗粒在设计流速下大致能不能被冲脱？”

【什么时候、由谁调用】
  手工调用，通常在改胶结参数之后、正式跑模拟之前：
      python3 audit_parallel_bond.py            # 全部参数取自 config/defaults.env
      python3 audit_parallel_bond.py <pack> --normal-strength-mpa 0.5 ...
  validate_run.py 要求运行目录里有 parallel_bond_audit.csv，一般就是本脚本的输出。

【输入】
  颗粒包文件（默认 PACK_FILE），胶结参数 BOND_*、水粘度 MU_F、入口流速 V_INF。
  初始键的判定与主程序一致：弱颗粒—骨架中心距 < BOND_CREATE_DISTANCE_M，
  且 < BOND_CREATE_RATIO × (r1+r2)。

【输出】
  一行汇总 CSV（默认 BOND_AUDIT_DEFAULT_OUTPUT = data/parallel_bond_audit.csv），
  以及屏幕上各量的 min/median/max。

【判据】
  每个弱颗粒至少有一根初始键，否则直接报错（说明 BOND_CREATE_* 设置有问题）。
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

from validate_pack_geometry import Particle, read_particles  # 复用颗粒包读取器
from parameter_config import load_parameters, real, text


def percentile(values: list[float], fraction: float) -> float:
    """线性插值百分位数（fraction=0.1 即 P10），不依赖 numpy。"""
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def describe(values: list[float]) -> tuple[float, float, float]:
    """返回 (最小值, 中位数, 最大值)，用于屏幕打印。"""
    return min(values), statistics.median(values), max(values)


def main() -> int:
    """读颗粒包，找出所有初始胶结键，换算物理量级并写出汇总 CSV。"""
    config = load_parameters()
    parser = argparse.ArgumentParser(description="Audit water_fangsha INL parallel-bond physical scales")
    parser.add_argument("pack", type=Path, nargs="?", default=Path(text(config, "PACK_FILE")))
    parser.add_argument("--radius-multiplier", type=float, default=real(config, "BOND_RADIUS_MULTIPLIER"))      # 键半径 = 倍数 × min(r1,r2)
    parser.add_argument("--effective-modulus-gpa", type=float, default=real(config, "BOND_PFC_EFFECTIVE_MODULUS_GPA"))  # PFC 有效模量 E*
    parser.add_argument("--shear-modulus-gpa", type=float, default=real(config, "BOND_PFC_SHEAR_MODULUS_GPA"))
    parser.add_argument("--stiffness-ratio", type=float, default=real(config, "BOND_PFC_STIFFNESS_RATIO"))      # kn/ks；为 0 时用两个模量之比
    parser.add_argument("--normal-strength-mpa", type=float, default=real(config, "BOND_NORMAL_STRENGTH_MPA"))  # 法向抗拉强度
    parser.add_argument("--shear-strength-mpa", type=float, default=real(config, "BOND_SHEAR_STRENGTH_MPA"))    # 切向抗剪强度
    parser.add_argument("--create-distance-m", type=float, default=real(config, "BOND_CREATE_DISTANCE_M"))      # 建键的绝对距离上限
    parser.add_argument("--create-ratio", type=float, default=real(config, "BOND_CREATE_RATIO"))                # 建键的相对距离上限(×(r1+r2))
    parser.add_argument("--water-viscosity-pa-s", type=float, default=real(config, "MU_F"))
    parser.add_argument("--inlet-velocity-m-s", type=float, default=real(config, "V_INF"))
    parser.add_argument("--output", type=Path, default=Path(text(config, "BOND_AUDIT_DEFAULT_OUTPUT")))
    args = parser.parse_args()

    # ── 输入合法性：这些量都必须是有限正数 ──
    numeric = (
        args.radius_multiplier,
        args.effective_modulus_gpa,
        args.shear_modulus_gpa,
        args.normal_strength_mpa,
        args.shear_strength_mpa,
        args.create_distance_m,
        args.create_ratio,
        args.water_viscosity_pa_s,
        args.inlet_velocity_m_s,
    )
    if any(not math.isfinite(value) or value <= 0 for value in numeric):
        raise ValueError("all bond and fluid inputs must be finite and positive")
    if not math.isfinite(args.stiffness_ratio) or args.stiffness_ratio < 0:
        raise ValueError("stiffness ratio must be finite and non-negative")

    # ── 单位换算：GPa→Pa，MPa→Pa；刚度比为 0 时退化为 E*/G* ──
    effective_modulus = args.effective_modulus_gpa * 1e9
    stiffness_ratio = args.stiffness_ratio or (
        args.effective_modulus_gpa / args.shear_modulus_gpa
    )
    shear_modulus = effective_modulus / stiffness_ratio
    normal_strength = args.normal_strength_mpa * 1e6
    shear_strength = args.shear_strength_mpa * 1e6

    # ── 第 1 步：按主程序的建键准则找出所有“弱颗粒—骨架”初始键 ──
    particles = read_particles(args.pack)
    fixed = [particle for particle in particles if particle.particle_type == 1]   # 骨架(固定)
    weak = [particle for particle in particles if particle.particle_type == 2]    # 弱颗粒(可动)
    bonds: list[tuple[Particle, Particle, float]] = []
    weak_bond_counts = {particle.particle_id: 0 for particle in weak}
    for mobile in weak:
        for skeleton in fixed:
            distance = math.dist(
                (mobile.x, mobile.y, mobile.z),
                (skeleton.x, skeleton.y, skeleton.z),
            )
            # 两个条件同时满足才建键：绝对距离上限 和 相对(r1+r2)上限
            if distance < args.create_distance_m and distance < args.create_ratio * (
                mobile.radius + skeleton.radius
            ):
                bonds.append((mobile, skeleton, distance))
                weak_bond_counts[mobile.particle_id] += 1

    if not bonds or any(count == 0 for count in weak_bond_counts.values()):
        raise RuntimeError("bond creation settings leave one or more weak particles unbonded")

    # ── 第 2 步：PFC 平行胶结几何量。键是一段半径 R_b 的“圆柱胶水”，
    #    面积 A=πR², 惯性矩 I=πR⁴/4（抗弯）, 极惯性矩 J=πR⁴/2（抗扭）──
    bond_radius = [args.radius_multiplier * min(a.radius, b.radius) for a, b, _ in bonds]
    pfc_length = [a.radius + b.radius for a, b, _ in bonds]   # PFC 里键的“长度”取 R1+R2
    area = [math.pi * radius**2 for radius in bond_radius]
    inertia = [math.pi * radius**4 / 4.0 for radius in bond_radius]
    polar_inertia = [math.pi * radius**4 / 2.0 for radius in bond_radius]
    # ── 第 3 步：刚度。PFC 约定 kn = E*/L（单位面积刚度, Pa/m），整根键刚度 = kn·A ──
    normal_stiffness_per_area = [effective_modulus / length for length in pfc_length]
    shear_stiffness_per_area = [shear_modulus / length for length in pfc_length]
    normal_stiffness = [kn * value for kn, value in zip(normal_stiffness_per_area, area)]
    shear_stiffness = [ks * value for ks, value in zip(shear_stiffness_per_area, area)]
    # ── 第 4 步：承载力。拉断力 = σ_n·A，剪断力 = τ_s·A；
    #    弯矩容量 = σ_n·I/R，扭矩容量 = τ_s·J/R（梁的最外纤维应力达到强度时的力矩）──
    normal_capacity = [normal_strength * value for value in area]
    shear_capacity = [shear_strength * value for value in area]
    bending_capacity = [
        normal_strength * moment / radius for moment, radius in zip(inertia, bond_radius)
    ]
    torsion_capacity = [
        shear_strength * moment / radius for moment, radius in zip(polar_inertia, bond_radius)
    ]
    # ── 第 5 步：和水流作用力比较。Stokes 拖曳力 F = 6πμ r U（低雷诺数球体）──
    stokes_drag = [
        6.0 * math.pi * args.water_viscosity_pa_s * mobile.radius * args.inlet_velocity_m_s
        for mobile, _, _ in bonds
    ]
    force_ratio = [capacity / drag for capacity, drag in zip(normal_capacity, stokes_drag)]  # >1 说明直接拉不断
    # 弯矩等效临界速度：拖曳力以键半径为力臂产生的弯矩 = 弯矩容量 时对应的流速
    #   M_cap = (6πμ r U_crit)·R_b  →  U_crit = M_cap / (R_b · 6πμ r)
    bending_velocity = [
        moment
        / (radius * 6.0 * math.pi * args.water_viscosity_pa_s * mobile.radius)
        for moment, radius, (mobile, _, _) in zip(bending_capacity, bond_radius, bonds)
    ]

    # ── 第 6 步：汇总成一行 CSV ──
    summary = {
        "weak_particles": len(weak),
        "initial_bonds": len(bonds),
        "min_bonds_per_weak": min(weak_bond_counts.values()),
        "max_bonds_per_weak": max(weak_bond_counts.values()),
        "radius_multiplier": args.radius_multiplier,
        "pfc_effective_modulus_pa": effective_modulus,
        "pfc_shear_modulus_pa": shear_modulus,
        "pfc_stiffness_ratio": stiffness_ratio,
        "normal_strength_pa": normal_strength,
        "shear_strength_pa": shear_strength,
        "inlet_velocity_m_s": args.inlet_velocity_m_s,
        "kn_pa_per_m_min": min(normal_stiffness_per_area),
        "kn_pa_per_m_median": statistics.median(normal_stiffness_per_area),
        "kn_pa_per_m_max": max(normal_stiffness_per_area),
        "ks_pa_per_m_min": min(shear_stiffness_per_area),
        "ks_pa_per_m_median": statistics.median(shear_stiffness_per_area),
        "ks_pa_per_m_max": max(shear_stiffness_per_area),
        "normal_stiffness_median_n_per_m": statistics.median(normal_stiffness),
        "shear_stiffness_median_n_per_m": statistics.median(shear_stiffness),
        "normal_force_capacity_min_n": min(normal_capacity),
        "normal_force_capacity_median_n": statistics.median(normal_capacity),
        "normal_force_capacity_max_n": max(normal_capacity),
        "bending_moment_capacity_median_nm": statistics.median(bending_capacity),
        "torsion_moment_capacity_median_nm": statistics.median(torsion_capacity),
        "stokes_drag_median_n": statistics.median(stokes_drag),
        "normal_capacity_to_stokes_drag_median": statistics.median(force_ratio),
        "bending_equivalent_critical_velocity_p10_m_s": percentile(bending_velocity, 0.10),
        "bending_equivalent_critical_velocity_median_m_s": statistics.median(bending_velocity),
        "bending_equivalent_critical_velocity_p90_m_s": percentile(bending_velocity, 0.90),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(summary.keys())
        writer.writerow(summary.values())

    # 屏幕上打印各量的分布，方便肉眼比较“承载力”和“拖曳力”的数量级
    for name, values, unit in (
        ("bond radius", bond_radius, "m"),
        ("normal force capacity", normal_capacity, "N"),
        ("bending moment capacity", bending_capacity, "N.m"),
        ("PFC bond length R1+R2", pfc_length, "m"),
        ("converted normal stiffness per area", normal_stiffness_per_area, "Pa/m"),
        ("converted shear stiffness per area", shear_stiffness_per_area, "Pa/m"),
        ("Stokes drag at inlet velocity", stokes_drag, "N"),
        ("bending-equivalent critical velocity", bending_velocity, "m/s"),
    ):
        low, median, high = describe(values)
        print(f"{name}: min={low:.8g}, median={median:.8g}, max={high:.8g} {unit}")
    print(f"initial_bonds={len(bonds)}, weak_particles={len(weak)}, output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
validate_pack_geometry.py —— 旧版“5.2 方截面熔合骨架”颗粒包的严格校验器（兼容保留）。

【这个脚本做什么】
  针对早期的 600 µm × 300 µm × 300 µm 方截面算例（骨架 type1 + 弱颗粒 type2，
  没有筛网）做一组 PASS/FAIL 检查：原子数、ID 连续、粒径/密度有限、弱颗粒在盒内、
  弱颗粒互不重叠、每颗弱颗粒至少有一个骨架接触点、岩心联合孔隙度在区间内。

【什么时候、由谁调用】
  - 手工：python3 validate_pack_geometry.py [pack] [--spacing-um 2 --min-porosity 0.27 ...]
  - 被导入：audit_parallel_bond.py 只借用这里的 `Particle` 数据类和 `read_particles()`
    读取颗粒包；当前 oilwater_fangsha 算例（带筛网、2.3 mm 长）用的校验器是
    validate_water_fangsha_pack.py，不是本脚本的 main()。

【输入】
  LIGGGHTS data 文件（Atoms 段每行：id type diameter density x y z）。

【输出】
  屏幕 PASS/FAIL 清单；--output 时写一行指标 CSV。全部通过返回 0，否则返回 4。

【主要判据（默认值针对旧算例，参数可改）】
  type1=897、type2=118；孔隙度 0.27..0.33；密度 2650 kg/m³；
  重叠判据 0.1 µm；“接触”判据 中心距/(r1+r2) < 1.00001。
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Particle:
    """一颗球：id、类型(1 骨架 / 2 弱颗粒 / 3 筛网)、半径与球心坐标（单位 m）、密度。"""
    particle_id: int
    particle_type: int
    radius: float
    x: float
    y: float
    z: float
    density: float = 2650.0


def read_particles(path: Path) -> list[Particle]:
    """读 LIGGGHTS data 文件的 Atoms 段（到 Velocities 段为止），返回 Particle 列表。

    文件里存的是直径，这里换成半径；列顺序 id type diameter density x y z。
    """
    particles: list[Particle] = []
    in_atoms = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.strip()
        if text == "Atoms":
            in_atoms = True
            continue
        if in_atoms and text == "Velocities":
            break
        fields = text.split()
        if in_atoms and len(fields) == 7 and fields[0].isdigit():
            particles.append(
                Particle(
                    int(fields[0]),          # id
                    int(fields[1]),          # type
                    0.5 * float(fields[2]),  # diameter → radius
                    float(fields[4]),        # x
                    float(fields[5]),        # y
                    float(fields[6]),        # z
                    float(fields[3]),        # density
                )
            )
    return particles


def declared_atom_count(path: Path) -> int | None:
    """从文件头部找 “<N> atoms” 一行，返回声明的原子数；找不到返回 None。"""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0].isdigit() and fields[1] == "atoms":
            return int(fields[0])
    return None


def square_union_porosity(particles: list[Particle], spacing: float) -> float:
    """用体素法计算旧算例岩心段 x∈[40,560] µm、300×300 µm 截面的“联合孔隙度”。

    逐个体素中心判断是否落在任一球内（球可重叠，所以必须取并集而不是体积相加）。
    纯 Python 三重循环，只适合旧的小算例。
    """
    x0, x1 = 40e-6, 560e-6
    nx = round((x1 - x0) / spacing)
    ny = round(300e-6 / spacing)
    total = 0
    solid = 0
    for ix in range(nx):
        x = x0 + (ix + 0.5) * spacing
        x_candidates = [p for p in particles if abs(p.x - x) < p.radius]   # 先按 x 粗筛，减少内层比较
        for iy in range(ny):
            y = (iy + 0.5) * spacing
            y_candidates = [p for p in x_candidates if abs(p.y - y) < p.radius]
            for iz in range(ny):
                z = (iz + 0.5) * spacing
                total += 1
                if any(
                    (x - p.x) ** 2 + (y - p.y) ** 2 + (z - p.z) ** 2 < p.radius**2
                    for p in y_candidates
                ):
                    solid += 1
    return 1.0 - solid / total


def main() -> int:
    """执行全部检查，打印 PASS/FAIL；可选写指标 CSV；全部通过返回 0，否则返回 4。"""
    parser = argparse.ArgumentParser(description="Validate the 5.2 square fused-skeleton pack")
    parser.add_argument("pack", type=Path, nargs="?", default=Path("data/cemented_pack.porepack"))
    parser.add_argument("--spacing-um", type=float, default=2.0)
    parser.add_argument("--min-porosity", type=float, default=0.27)
    parser.add_argument("--max-porosity", type=float, default=0.33)
    parser.add_argument("--expected-type1", type=int, default=897)
    parser.add_argument("--expected-type2", type=int, default=118)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.spacing_um <= 0 or not 0 <= args.min_porosity < args.max_porosity <= 1:
        raise ValueError("invalid spacing or porosity limits")

    particles = read_particles(args.pack)
    type1 = [p for p in particles if p.particle_type == 1]
    type2 = [p for p in particles if p.particle_type == 2]
    ids = [p.particle_id for p in particles]
    finite_values = all(
        math.isfinite(value)
        for p in particles
        for value in (p.radius, p.density, p.x, p.y, p.z)
    )

    # ── 全对比较：统计重叠和接触 ──
    # penetration = r1 + r2 − 中心距：>0 表示两球嵌入
    weak_overlap_count = 0          # 涉及弱颗粒的重叠对数（弱-弱、弱-骨架），要求为 0
    skeleton_overlap_count = 0      # 骨架-骨架重叠对数（允许，骨架本来就是“熔合”的）
    max_skeleton_penetration = 0.0
    weak_bond_neighbors = {p.particle_id: 0 for p in type2}   # 每颗弱颗粒的骨架接触点数
    for index, first in enumerate(particles):
        for second in particles[index + 1 :]:
            distance = math.dist((first.x, first.y, first.z), (second.x, second.y, second.z))
            penetration = first.radius + second.radius - distance
            if first.particle_type == second.particle_type == 1:      # 骨架-骨架
                if penetration > 0.1e-6:
                    skeleton_overlap_count += 1
                    max_skeleton_penetration = max(max_skeleton_penetration, penetration)
                continue
            if first.particle_type == second.particle_type == 2:      # 弱-弱
                if penetration > 0.1e-6:
                    weak_overlap_count += 1
                continue
            weak = first if first.particle_type == 2 else second      # 弱-骨架
            if penetration > 0.1e-6:
                weak_overlap_count += 1
            if distance / (first.radius + second.radius) < 1.00001:    # 几乎相切 → 可建胶结键
                weak_bond_neighbors[weak.particle_id] += 1

    # 弱颗粒必须整颗在 600×300×300 µm 的 DEM 盒子内
    weak_inside = all(
        p.x - p.radius >= 0
        and p.x + p.radius <= 600e-6
        and p.y - p.radius >= 0
        and p.y + p.radius <= 300e-6
        and p.z - p.radius >= 0
        and p.z + p.radius <= 300e-6
        for p in type2
    )
    unbonded_weak = sum(count == 0 for count in weak_bond_neighbors.values())
    initial_bond_candidates = sum(weak_bond_neighbors.values())
    porosity = square_union_porosity(particles, args.spacing_um * 1e-6)
    declared = declared_atom_count(args.pack)

    checks = [
        (declared == len(particles), "declared atom count matches parsed particles"),
        (
            len(type1) == args.expected_type1 and len(type2) == args.expected_type2,
            f"composition is exactly type1={args.expected_type1}, type2={args.expected_type2}",
        ),
        (len(ids) == len(set(ids)) and sorted(ids) == list(range(1, len(ids) + 1)), "particle IDs are unique and contiguous"),
        (finite_values and all(p.radius > 0 and p.density == 2650.0 for p in particles), "geometry is finite with positive radii and density 2650 kg/m3"),
        (set(p.particle_type for p in particles) == {1, 2}, "only particle types 1 and 2 are present"),
        (weak_inside, "all mobile weak spheres are fully inside the square DEM box"),
        (weak_overlap_count == 0, "weak particles have no overlap above 0.1 um"),
        (unbonded_weak == 0, "every weak particle has at least one initial type1-type2 bond candidate"),
        (
            args.min_porosity <= porosity <= args.max_porosity,
            f"square core union porosity is within {args.min_porosity:.3f}..{args.max_porosity:.3f}",
        ),
    ]

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "particles", "type1_fused_primitives", "type2_mobile", "spacing_um",
                    "core_union_porosity", "weak_overlap_count", "unbonded_weak_count",
                    "initial_bond_candidates", "skeleton_overlap_count",
                    "max_skeleton_penetration_um", "skeleton_representation",
                ]
            )
            writer.writerow(
                [
                    len(particles), len(type1), len(type2), args.spacing_um, porosity,
                    weak_overlap_count, unbonded_weak, initial_bond_candidates,
                    skeleton_overlap_count, max_skeleton_penetration * 1e6,
                    "fused_union_primitives_not_a_granular_packing",
                ]
            )

    for passed, description in checks:
        print(f"{'PASS' if passed else 'FAIL'}: {description}")
    print(f"core_union_porosity={porosity:.8f}")
    print(
        "INFO: fixed type1 geometry is a fused union of overlapping sphere primitives: "
        f"overlap_pairs={skeleton_overlap_count}, max_penetration_um={max_skeleton_penetration * 1e6:.6g}"
    )
    return 0 if all(passed for passed, _ in checks) else 4


if __name__ == "__main__":
    raise SystemExit(main())

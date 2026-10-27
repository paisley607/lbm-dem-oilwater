#!/usr/bin/env python3
"""
validate_water_fangsha_pack.py —— 离线校验 oilwater_fangsha 颗粒包：骨架、细颗粒、出口筛网。

【这个脚本做什么】
  在真正跑 LBM-DEM 之前，用纯几何方法检查 make_water_fangsha_pack.py 生成的颗粒包
  是否“像一块能出砂的岩心 + 一张合理的筛网”：
    - 计算域尺寸与配置一致，只含 type 1/2/3
    - 骨架粒径在允许范围、细颗粒直径等于设定值
    - 骨架联合孔隙度（体素法）落在 PHI_MIN..PHI_MAX
    - 细颗粒浓度（体积分数）落在 FINE_CONCENTRATION_MIN..MAX
    - ≥99% 的细颗粒与某颗骨架球相切（否则主程序建不了胶结键）
    - 骨架最大重叠比例不超过 SKELETON_MAX_OVERLAP_FRACTION
    - 孔隙从入口到出口连通（六邻域连通分量）
    - 筛网开孔率在 (0,1) 内，且“最大通过孔径 / 细颗粒直径”落在
      SCREEN_MIN_APERTURE_RATIO..SCREEN_MAX_APERTURE_RATIO
  任一条不满足就抛 ValueError（退出码非 0）。

【什么时候、由谁调用】
  - 自动：run_production.sh 每次启动时先调用本脚本，失败则不编译、不运行主程序。
  - 手工：python3 validate_water_fangsha_pack.py data/oilwater_fangsha.porepack --output xxx.csv
    （几何参数可用命令行覆盖；--spacing-um 应等于 LBM 格距 = REFERENCE_DIAMETER_UM / N）

【输入】
  .porepack 文件（Atoms 段每行 id type diameter density x y z）+ config/defaults.env。

【输出】
  屏幕打印全部指标；--output 时写一行指标 CSV（run_production.sh 写到 $OUTDIR/pack_geometry.csv）；
  最后一行 "[PASS] water_fangsha rock, fines, and screen geometry"。
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy import ndimage
from parameter_config import load_parameters, real


def read_pack(path):
    """读 .porepack：返回 (三个方向的盒子边界, 颗粒列表)。

    每个颗粒是元组 (id, type, radius, density, x, y, z)；文件里存直径，这里换成半径。
    """
    lines = path.read_text(encoding="ascii").splitlines()
    atom_line = next(line for line in lines if line.endswith(" atoms"))
    n = int(atom_line.split()[0])                         # 头部声明的原子数
    bounds = []
    for axis in ("x", "y", "z"):                          # 找 "lo hi xlo xhi" 这样的边界行
        line = next(line for line in lines if len(line.split()) >= 4 and line.split()[2] == f"{axis}lo" and line.split()[3] == f"{axis}hi")
        values = line.split()
        bounds.append((float(values[0]), float(values[1])))
    start = lines.index("Atoms") + 2                      # "Atoms" 行 + 一个空行之后是数据
    particles = []
    for line in lines[start:]:
        if not line.strip():
            continue
        values = line.split()
        if len(values) < 7:
            continue
        particles.append((int(values[0]), int(values[1]), 0.5 * float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6])))
    if len(particles) != n:
        raise ValueError(f"header says {n} atoms but read {len(particles)}")
    return bounds, particles


def union_mask(particles, spacing, x0, x1, ly, lz):
    """把所有球体素化到岩心段 [x0,x1]×[0,ly]×[0,lz]，返回 (固体掩膜, 体素总数)。

    体素分辨率 spacing 取 LBM 格距，这样算出的孔隙度和主程序 IBM 的口径一致。
    """
    nx, ny, nz = round((x1 - x0) / spacing), round(ly / spacing), round(lz / spacing)
    mask = np.zeros((nx, ny, nz), dtype=bool)
    for _, _, r, _, px, py, pz in particles:
        # 球的包围盒 → 体素索引范围（裁剪到掩膜内）
        ix0 = max(0, int(math.floor((px - r - x0) / spacing)))
        ix1 = min(nx - 1, int(math.floor((px + r - x0) / spacing)))
        iy0 = max(0, int(math.floor((py - r) / spacing)))
        iy1 = min(ny - 1, int(math.floor((py + r) / spacing)))
        iz0 = max(0, int(math.floor((pz - r) / spacing)))
        iz1 = min(nz - 1, int(math.floor((pz + r) / spacing)))
        xs = x0 + (np.arange(ix0, ix1 + 1) + 0.5) * spacing   # 体素中心坐标
        ys = (np.arange(iy0, iy1 + 1) + 0.5) * spacing
        zs = (np.arange(iz0, iz1 + 1) + 0.5) * spacing
        inside = (xs[:, None, None] - px) ** 2 + (ys[None, :, None] - py) ** 2 + (zs[None, None, :] - pz) ** 2 <= r * r
        mask[ix0 : ix1 + 1, iy0 : iy1 + 1, iz0 : iz1 + 1] |= inside   # 并集：重叠部分只算一次
    return mask, mask.size


def screen_mask(particles, spacing, ly, lz, screen_x):
    """把筛网珠子投影到 x=screen_x 的 y-z 平面，返回二维固体掩膜（True=被丝挡住）。"""
    ny, nz = round(ly / spacing), round(lz / spacing)
    mask = np.zeros((ny, nz), dtype=bool)
    ys = (np.arange(ny) + 0.5) * spacing
    zs = (np.arange(nz) + 0.5) * spacing
    for _, _, radius, _, px, py, pz in particles:
        if abs(px - screen_x) > 1e-10:                    # 所有筛网珠子必须严格在同一平面上
            raise ValueError(f"screen bead is not on x={screen_x}: x={px}")
        mask |= (ys[:, None] - py) ** 2 + (zs[None, :] - pz) ** 2 <= radius**2
    return mask


def main():
    """读颗粒包 → 算各项几何指标 → 逐项判据 → 打印/写 CSV。"""
    config = load_parameters()
    parser = argparse.ArgumentParser()
    parser.add_argument("pack", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--min-porosity", type=float, default=real(config, "PHI_MIN"))
    parser.add_argument("--max-porosity", type=float, default=real(config, "PHI_MAX"))
    parser.add_argument("--spacing-um", type=float, default=real(config, "REFERENCE_DIAMETER_UM") / real(config, "N"))   # 岩心体素 = LBM 格距
    parser.add_argument("--screen-spacing-um", type=float, default=real(config, "SCREEN_VALIDATION_SPACING_UM"))         # 筛网平面用更细的像素
    parser.add_argument("--lx-mm", type=float, default=real(config, "DOMAIN_LX_M") * 1000.0)
    parser.add_argument("--ly-mm", type=float, default=real(config, "DOMAIN_LY_M") * 1000.0)
    parser.add_argument("--lz-mm", type=float, default=real(config, "DOMAIN_LZ_M") * 1000.0)
    parser.add_argument("--core-x0-mm", type=float, default=real(config, "CORE_X0_M") * 1000.0)
    parser.add_argument("--core-x1-mm", type=float, default=real(config, "CORE_X1_M") * 1000.0)
    parser.add_argument("--max-skeleton-overlap", type=float, default=real(config, "SKELETON_MAX_OVERLAP_FRACTION"))
    parser.add_argument("--screen-x-mm", type=float, default=real(config, "SCREEN_X_M") * 1000.0)
    parser.add_argument("--screen-wire-diameter-um", type=float, default=real(config, "SCREEN_WIRE_DIAMETER_UM"))
    parser.add_argument("--fine-diameter-um", type=float, default=real(config, "FINE_DIAMETER_UM"))
    parser.add_argument("--min-screen-aperture-ratio", type=float, default=real(config, "SCREEN_MIN_APERTURE_RATIO"))
    parser.add_argument("--max-screen-aperture-ratio", type=float, default=real(config, "SCREEN_MAX_APERTURE_RATIO"))
    args = parser.parse_args()

    # ── 第 1 步：读文件，核对计算域尺寸和颗粒类型 ──
    bounds, particles = read_pack(args.pack)
    expected = ((0.0, args.lx_mm * 1e-3), (0.0, args.ly_mm * 1e-3), (0.0, args.lz_mm * 1e-3))
    if any(abs(a - b) > 1e-12 or abs(c - d) > 1e-12 for (a, c), (b, d) in zip(bounds, expected)):
        raise ValueError(f"domain mismatch: {bounds}")
    skeleton = [p for p in particles if p[1] == 1]
    fines = [p for p in particles if p[1] == 2]
    screen = [p for p in particles if p[1] == 3]
    if not skeleton or not fines or not screen:
        raise ValueError("skeleton, fine particles, and type-3 screen beads are required")
    if any(p[1] not in (1, 2, 3) for p in particles):
        raise ValueError("only particle types 1=skeleton, 2=fines, 3=screen are allowed")

    # ── 第 2 步：粒径检查 ──
    diameters = [2.0 * p[2] for p in skeleton]
    fine_diameters = [2.0 * p[2] for p in fines]
    skeleton_min = real(config, "SKELETON_DIAMETER_MIN_MM") * 1e-3
    skeleton_max = real(config, "SKELETON_DIAMETER_MAX_MM") * 1e-3
    if min(diameters) < skeleton_min - 1e-9 or max(diameters) > skeleton_max + 1e-9:
        raise ValueError(f"skeleton diameter is outside {skeleton_min * 1e3:g}--{skeleton_max * 1e3:g} mm")
    expected_fine_diameter = args.fine_diameter_um * 1e-6
    if max(abs(d - expected_fine_diameter) for d in fine_diameters) > 1e-9:
        raise ValueError(f"fine diameter is not {args.fine_diameter_um:g} um")

    # ── 第 3 步：骨架联合孔隙度与细颗粒浓度（体素法）──
    spacing = args.spacing_um * 1e-6
    mask, total = union_mask(skeleton, spacing, args.core_x0_mm * 1e-3, args.core_x1_mm * 1e-3, args.ly_mm * 1e-3, args.lz_mm * 1e-3)
    solid_voxels = int(np.count_nonzero(mask))
    skeleton_porosity = 1.0 - solid_voxels / total
    fine_volume = len(fines) * (4.0 * math.pi / 3.0) * (0.5 * expected_fine_diameter) ** 3
    skeleton_union_volume = solid_voxels * spacing**3
    fine_concentration = fine_volume / (skeleton_union_volume + fine_volume)   # 细颗粒占固相的体积分数

    # ── 第 4 步：筛网孔径计算 ──
    # 筛网由相切珠子排成的丝组成，真实的“孔”不是简单的矩形，所以不用名义孔径，
    # 而是在 y-z 平面上做二维投影后用距离变换找“最大能通过的圆”：
    #   (a) screen_mask 把珠子投影成像素（True=丝）；
    #   (b) distance_transform_edt 对每个开孔像素算它到最近丝像素的距离；
    #   (c) 最大距离 = 孔中能放下的最大内切圆半径，乘 2 得最大通过孔径；
    #   (d) 孔径 / 细颗粒直径 = aperture_ratio：<1 细颗粒根本过不去，太大则筛网形同虚设，
    #       所以要求落在 SCREEN_MIN_APERTURE_RATIO..SCREEN_MAX_APERTURE_RATIO。
    screen_x = args.screen_x_mm * 1e-3
    expected_screen_radius = 0.5 * args.screen_wire_diameter_um * 1e-6
    if any(abs(p[2] - expected_screen_radius) > 1e-10 for p in screen):
        raise ValueError("screen bead diameter does not match --screen-wire-diameter-um")
    screen_spacing = args.screen_spacing_um * 1e-6
    if screen_spacing <= 0.0:
        raise ValueError("screen spacing must be positive")
    smask = screen_mask(screen, screen_spacing, args.ly_mm * 1e-3, args.lz_mm * 1e-3, screen_x)
    screen_open_fraction = 1.0 - float(np.mean(smask))                     # 开孔率
    clearance = ndimage.distance_transform_edt(~smask) * screen_spacing     # 每个开孔像素到最近丝的距离 [m]
    maximum_opening = 2.0 * float(np.max(clearance))                        # 最大内切圆直径 = 最大通过孔径
    aperture_ratio = maximum_opening / min(fine_diameters)

    # ── 第 5 步：孔隙连通性。把孔隙体素做六邻域连通标记，看入口面和出口面是否共享同一标签 ──
    pore_labels, _ = ndimage.label(~mask, structure=ndimage.generate_binary_structure(3, 1))
    inlet_labels = np.unique(pore_labels[0, :, :])
    outlet_labels = np.unique(pore_labels[-1, :, :])
    through_labels = np.intersect1d(inlet_labels[inlet_labels > 0], outlet_labels[outlet_labels > 0])

    # ── 第 6 步：细颗粒—骨架接触。中心距/(r1+r2) 与 1 的偏差 <= BOND_CONTACT_TOLERANCE 视为相切 ──
    bonded_fines = 0
    bond_contact_tolerance = real(config, "BOND_CONTACT_TOLERANCE")
    for fine in fines:
        if any(abs(math.dist((fine[4], fine[5], fine[6]), (skel[4], skel[5], skel[6])) / (fine[2] + skel[2]) - 1.0) <= bond_contact_tolerance for skel in skeleton):
            bonded_fines += 1

    # ── 第 7 步：骨架之间的重叠比例 = 1 − 中心距/(r1+r2)，统计对数和最大值 ──
    overlap_pairs = 0
    max_overlap = 0.0
    for index, first in enumerate(skeleton):
        for second in skeleton[index + 1 :]:
            separation = math.dist((first[4], first[5], first[6]), (second[4], second[5], second[6]))
            overlap = max(0.0, 1.0 - separation / (first[2] + second[2]))
            if overlap > 0.0:
                overlap_pairs += 1
                max_overlap = max(max_overlap, overlap)

    metrics = {
        "atoms": len(particles), "skeleton": len(skeleton), "fines": len(fines), "screen": len(screen),
        "skeleton_union_porosity": skeleton_porosity,
        "fine_concentration_vs_skeleton_union": fine_concentration,
        "fines_with_contact_candidate": bonded_fines,
        "skeleton_diameter_min_mm": min(diameters) * 1e3,
        "skeleton_diameter_max_mm": max(diameters) * 1e3,
        "skeleton_overlap_pairs": overlap_pairs,
        "skeleton_max_overlap_fraction": max_overlap,
        "through_pore_components": len(through_labels),
        "screen_x_mm": screen_x * 1e3,
        "screen_wire_diameter_um": 2.0 * expected_screen_radius * 1e6,
        "screen_open_fraction": screen_open_fraction,
        "screen_maximum_opening_um": maximum_opening * 1e6,
        "screen_aperture_to_fine_ratio": aperture_ratio,
    }

    # ── 第 8 步：逐项判据，任一不满足即抛错 ──
    if not args.min_porosity <= skeleton_porosity <= args.max_porosity:
        raise ValueError(f"skeleton union porosity {skeleton_porosity:.6f} outside {args.min_porosity}..{args.max_porosity}")
    if bonded_fines < int(0.99 * len(fines)):
        raise ValueError(f"only {bonded_fines}/{len(fines)} fines have a skeleton contact")
    if max_overlap > args.max_skeleton_overlap + 1e-8:
        raise ValueError(f"maximum skeleton overlap {max_overlap:.6f} exceeds {args.max_skeleton_overlap:.6f}")
    if not len(through_labels):
        raise ValueError("no six-neighbor pore component crosses the core from inlet to outlet")
    fine_min = real(config, "FINE_CONCENTRATION_MIN")
    fine_max = real(config, "FINE_CONCENTRATION_MAX")
    if not fine_min <= fine_concentration <= fine_max:
        raise ValueError(f"fine concentration {fine_concentration:.6f} outside {fine_min:g}..{fine_max:g}")
    if not 0.0 < screen_open_fraction < 1.0:
        raise ValueError(f"invalid screen open fraction {screen_open_fraction:.6f}")
    if not args.min_screen_aperture_ratio <= aperture_ratio <= args.max_screen_aperture_ratio:
        raise ValueError(
            f"screen aperture/fine ratio {aperture_ratio:.6f} outside "
            f"{args.min_screen_aperture_ratio}..{args.max_screen_aperture_ratio}"
        )

    # ── 输出 ──
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="ascii") as handle:
            writer = csv.writer(handle)
            writer.writerow(metrics.keys())
            writer.writerow(metrics.values())
    for key, value in metrics.items():
        print(f"{key}={value}")
    print("[PASS] water_fangsha rock, fines, and screen geometry")


if __name__ == "__main__":
    main()

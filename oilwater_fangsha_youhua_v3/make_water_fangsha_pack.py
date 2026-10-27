#!/usr/bin/env python3
"""
make_water_fangsha_pack.py —— 生成出砂模拟用的“骨架 + 细颗粒 + 出口筛网”颗粒包。

【这个脚本做什么】
  按 config/defaults.env 里的几何与粒径参数，生成一个 LIGGGHTS 可读的 .porepack
  数据文件，里面有三类球：
    type 1  骨架颗粒（固定不动，相互允许少量重叠，合起来当作“胶结岩石”）
    type 2  细颗粒 / 弱颗粒（每颗贴在一颗骨架颗粒表面，运行时靠平行胶结键固定，
            是唯一可能被水流冲走、产出的颗粒）
    type 3  筛网珠子（固定不动，排成正交的两组“丝”，模拟出口防砂筛网）

【什么时候、由谁调用】
  手工调用，只在需要重新生成或改变颗粒包时执行：
      python3 make_water_fangsha_pack.py data/oilwater_fangsha.porepack [--seed 1 ...]
  生成后应立即用 validate_water_fangsha_pack.py 校验；run_production.sh 每次启动
  前也会再校验一遍。

【输入】
  全部参数默认来自 config/defaults.env（DOMAIN_L*_M、CORE_X*_M、REFERENCE_CORE_POROSITY、
  FINE_FRACTION、SKELETON_DIAMETER_*、SCREEN_* 等），命令行可覆盖。

【输出】
  一个 LIGGGHTS data 文件（Atoms 段每行：id type diameter density x y z），
  屏幕打印骨架数、联合孔隙度、细颗粒数与浓度、筛网珠子数、名义孔径等。

【生成思路（详见各函数）】
  1. build_skeleton：先按 FCC 密排格点摆骨架球（加小抖动和粒径分布），再随机补洞，
     直到体素联合孔隙度降到目标值。
  2. build_fines：每颗细颗粒放在随机骨架球表面的随机方向上，与该骨架球恰好相切，
     且不与其他骨架球/细颗粒重叠；数量由“细颗粒体积分数”反推。
  3. build_screen：在 x=screen_x 平面上，用相切的珠子排出 y、z 两个方向的丝，
     得到方孔筛网。
"""

from __future__ import annotations

import argparse
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from parameter_config import load_parameters, real, text


class SpatialBins:
    """简易空间哈希（格子分桶）：查“附近有哪些球”时只看相邻格子，避免 O(N²) 的全对比较。"""

    def __init__(self, cell_size: float) -> None:
        self.cell_size = cell_size
        self.bins: dict[tuple[int, int, int], list[int]] = defaultdict(list)

    def key(self, x: float, y: float, z: float) -> tuple[int, int, int]:
        """把坐标换成所在格子的整数索引。"""
        return (int(math.floor(x / self.cell_size)), int(math.floor(y / self.cell_size)), int(math.floor(z / self.cell_size)))

    def add(self, index: int, x: float, y: float, z: float) -> None:
        """把第 index 个球登记到它所在的格子。"""
        self.bins[self.key(x, y, z)].append(index)

    def nearby(self, x: float, y: float, z: float, reach: float) -> list[int]:
        """返回距离 (x,y,z) 在 reach 以内可能相交的所有球的索引（取相邻若干层格子）。"""
        cx, cy, cz = self.key(x, y, z)
        width = int(math.ceil(reach / self.cell_size))
        result: list[int] = []
        for ix in range(cx - width, cx + width + 1):
            for iy in range(cy - width, cy + width + 1):
                for iz in range(cz - width, cz + width + 1):
                    result.extend(self.bins.get((ix, iy, iz), ()))
        return result


def draw_diameter(rng: random.Random, minimum: float, maximum: float, mean: float, sigma: float) -> float:
    """从截断正态分布抽一个骨架粒径：超出 [minimum, maximum] 就重抽。"""
    while True:
        value = rng.gauss(mean, sigma)
        if minimum <= value <= maximum:
            return value


def add_sphere_to_mask(mask: np.ndarray, spacing: float, x0: float, x: float, y: float, z: float, radius: float) -> int:
    """把一个球“并”进体素掩膜（True=固体），返回新增的固体体素数。

    体素掩膜只覆盖岩心段 [x0, x1]，分辨率 spacing 与 LBM 格子一致，
    这样算出的“联合孔隙度”和主程序 IBM 看到的孔隙度口径相同。
    """
    nx, ny, nz = mask.shape
    # 球的包围盒在体素索引上的范围（裁剪到掩膜边界内）
    ix0 = max(0, int(math.floor((x - radius - x0) / spacing)))
    ix1 = min(nx - 1, int(math.floor((x + radius - x0) / spacing)))
    iy0 = max(0, int(math.floor((y - radius) / spacing)))
    iy1 = min(ny - 1, int(math.floor((y + radius) / spacing)))
    iz0 = max(0, int(math.floor((z - radius) / spacing)))
    iz1 = min(nz - 1, int(math.floor((z + radius) / spacing)))
    # 体素中心坐标
    xs = x0 + (np.arange(ix0, ix1 + 1) + 0.5) * spacing
    ys = (np.arange(iy0, iy1 + 1) + 0.5) * spacing
    zs = (np.arange(iz0, iz1 + 1) + 0.5) * spacing
    # 体素中心到球心距离 <= 半径 → 该体素属于球
    inside = (xs[:, None, None] - x) ** 2 + (ys[None, :, None] - y) ** 2 + (zs[None, None, :] - z) ** 2 <= radius**2
    region = mask[ix0 : ix1 + 1, iy0 : iy1 + 1, iz0 : iz1 + 1]
    new_voxels = int(np.count_nonzero(inside & ~region))   # 原来是孔隙、现在变固体的体素
    region |= inside
    return new_voxels


def overlaps(candidate: tuple[float, float, float, float], particles: list[tuple[int, int, float, float, float, float, float]], bins: SpatialBins, max_overlap: float, maximum_diameter: float) -> bool:
    """判断候选球是否与已有球重叠超过允许比例。

    允许重叠：中心距 >= (r1+r2)·(1-max_overlap)。骨架颗粒允许少量互相嵌入，
    代表颗粒之间被胶结物“焊”在一起，而不是松散堆积。
    """
    radius, x, y, z = candidate
    for index in bins.nearby(x, y, z, radius + maximum_diameter):
        other = particles[index]   # 元组 (id, type, radius, x, y, z, density)
        distance = math.dist((x, y, z), (other[3], other[4], other[5]))
        if distance < (radius + other[2]) * (1.0 - max_overlap):
            return True
    return False


def build_skeleton(rng: random.Random, core_x0: float, core_x1: float, ly: float, lz: float, spacing: float, target_porosity: float, max_overlap: float, lattice_a: float, diameter_min: float, diameter_max: float, diameter_mean: float, diameter_sigma: float, density: float) -> tuple[list[tuple[int, int, float, float, float, float, float]], np.ndarray]:
    """在岩心段 [core_x0, core_x1]×[0,ly]×[0,lz] 内生成固定骨架球，直到联合孔隙度 <= 目标值。

    返回 (骨架球列表, 固体体素掩膜)。每个球是元组 (id, type=1, radius, x, y, z, density)。
    """
    # 体素掩膜：用来实时监控“联合孔隙度”（球可以重叠，所以不能简单把球体积相加）
    nx = round((core_x1 - core_x0) / spacing)
    ny = round(ly / spacing)
    nz = round(lz / spacing)
    mask = np.zeros((nx, ny, nz), dtype=bool)
    skeleton: list[tuple[int, int, float, float, float, float, float]] = []
    bins = SpatialBins(diameter_max)
    target_solid = 1.0 - target_porosity        # 目标固体体积分数

    # ── 阶段 1：FCC 密排格点作为起始骨架 ──
    # 纯随机顺序投放(RSA)大约在孔隙度 50% 就“堵死”了，到不了 0.35~0.38 的目标。
    # 所以先用面心立方(FCC)格点铺一层接近相切的致密骨架：
    # 每个晶胞 4 个基点，晶格常数 lattice_a；再加小抖动和粒径分布打破规则感。
    basis = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    candidates: list[tuple[float, float, float]] = []
    ix_min = int(math.floor(core_x0 / lattice_a)) - 1
    ix_max = int(math.ceil(core_x1 / lattice_a)) + 1
    iy_max = int(math.ceil(ly / lattice_a)) + 1
    iz_max = int(math.ceil(lz / lattice_a)) + 1
    for ix in range(ix_min, ix_max + 1):
        for iy in range(-1, iy_max + 1):
            for iz in range(-1, iz_max + 1):
                for bx, by, bz in basis:
                    candidates.append(((ix + bx) * lattice_a, (iy + by) * lattice_a, (iz + bz) * lattice_a))
    rng.shuffle(candidates)   # 随机顺序投放，避免“从一头长过去”的方向性

    def try_add(x: float, y: float, z: float, allowed_overlap: float) -> bool:
        """在 (x,y,z) 试放一颗随机粒径的骨架球：出界或重叠超限则放弃。"""
        diameter = draw_diameter(rng, diameter_min, diameter_max, diameter_mean, diameter_sigma)
        radius = 0.5 * diameter
        if not (core_x0 + radius <= x <= core_x1 - radius and radius <= y <= ly - radius and radius <= z <= lz - radius):
            return False   # 整颗球必须完全落在岩心段内
        if overlaps((radius, x, y, z), skeleton, bins, allowed_overlap, diameter_max):
            return False
        add_sphere_to_mask(mask, spacing, core_x0, x, y, z, radius)
        skeleton.append((len(skeleton) + 1, 1, radius, x, y, z, density))
        bins.add(len(skeleton) - 1, x, y, z)
        return True

    for x, y, z in candidates:
        jitter = 0.008e-3   # 格点抖动 ±8 µm
        try_add(x + rng.uniform(-jitter, jitter), y + rng.uniform(-jitter, jitter), z + rng.uniform(-jitter, jitter), max_overlap)
        if float(np.mean(mask)) >= target_solid:   # 固体分数到达目标 → 孔隙度达标，提前结束
            return skeleton, mask

    # ── 阶段 2：格点铺完仍未达标，就随机补洞（填充局部缺陷）──
    for _ in range(250000):
        diameter = draw_diameter(rng, diameter_min, diameter_max, diameter_mean, diameter_sigma)
        radius = 0.5 * diameter
        x = rng.uniform(core_x0 + radius, core_x1 - radius)
        y = rng.uniform(radius, ly - radius)
        z = rng.uniform(radius, lz - radius)
        if overlaps((radius, x, y, z), skeleton, bins, max_overlap, diameter_max):
            continue
        add_sphere_to_mask(mask, spacing, core_x0, x, y, z, radius)
        skeleton.append((len(skeleton) + 1, 1, radius, x, y, z, density))
        bins.add(len(skeleton) - 1, x, y, z)
        if float(np.mean(mask)) >= target_solid:
            return skeleton, mask
    achieved_porosity = 1.0 - float(np.mean(mask))
    raise RuntimeError(
        f"could not reach target porosity {target_porosity:.3f}; achieved {achieved_porosity:.3f} with controlled overlap <= {max_overlap:.3f}"
    )


def build_fines(rng: random.Random, skeleton: list[tuple[int, int, float, float, float, float, float]], skeleton_mask: np.ndarray, spacing: float, core_x0: float, lx: float, ly: float, lz: float, fine_fraction: float, fine_diameter: float, density: float, skeleton_diameter_max: float) -> list[tuple[int, int, float, float, float, float, float]]:
    """把细颗粒(type 2)一颗颗贴到骨架球表面，直到达到目标细颗粒体积分数。

    每颗细颗粒与“母骨架球”恰好相切（这样主程序一启动就能按距离准则建立胶结键），
    同时不能嵌进其他骨架球，也不能与已放的细颗粒重叠。
    """
    fine_radius = 0.5 * fine_diameter
    # 细颗粒个数：由 fine_fraction = V_fine / (V_fine + V_skeleton) 反推
    skeleton_union_volume = float(np.count_nonzero(skeleton_mask)) * spacing**3
    fine_volume = 4.0 * math.pi * fine_radius**3 / 3.0
    target_fine_volume = fine_fraction / (1.0 - fine_fraction) * skeleton_union_volume
    fine_count = max(1, int(round(target_fine_volume / fine_volume)))
    # 两个空间哈希：一个查骨架、一个查已放的细颗粒
    skeleton_bins = SpatialBins(max(particle[2] for particle in skeleton) * 2.0)
    for index, particle in enumerate(skeleton):
        skeleton_bins.add(index, particle[3], particle[4], particle[5])
    fine_bins = SpatialBins(2.0 * fine_radius)
    fines: list[tuple[int, int, float, float, float, float, float]] = []

    for _ in range(fine_count):
        for _attempt in range(6000):   # 每颗最多试 6000 个位置
            # (1) 随机挑一颗母骨架球，在球面上随机取一个方向（均匀分布于球面）
            parent_index = rng.randrange(len(skeleton))
            parent = skeleton[parent_index]
            theta = rng.uniform(0.0, 2.0 * math.pi)
            axial = rng.uniform(-1.0, 1.0)
            radial = math.sqrt(max(0.0, 1.0 - axial * axial))
            direction = (radial * math.cos(theta), radial * math.sin(theta), axial)
            # (2) 沿该方向放到“恰好相切”的位置（乘 1+1e-8 避免数值上微小嵌入）
            distance = (parent[2] + fine_radius) * (1.0 + 1.0e-8)
            x = parent[3] + distance * direction[0]
            y = parent[4] + distance * direction[1]
            z = parent[5] + distance * direction[2]
            # (3) 细颗粒必须整个在计算域内（注意用全域 lx，不是岩心段）
            if not (fine_radius <= x <= lx - fine_radius and fine_radius <= y <= ly - fine_radius and fine_radius <= z <= lz - fine_radius):
                continue
            # (4) 不能嵌进母球之外的其他骨架球
            blocked = False
            for index in skeleton_bins.nearby(x, y, z, skeleton_diameter_max):
                other = skeleton[index]
                if index != parent_index and math.dist((x, y, z), (other[3], other[4], other[5])) < other[2] + fine_radius:
                    blocked = True
                    break
            if blocked:
                continue
            # (5) 不能与已放的细颗粒重叠
            for index in fine_bins.nearby(x, y, z, 2.0 * fine_radius):
                other = fines[index]
                if math.dist((x, y, z), (other[3], other[4], other[5])) < 2.0 * fine_radius:
                    blocked = True
                    break
            if blocked:
                continue
            fine = (len(skeleton) + len(fines) + 1, 2, fine_radius, x, y, z, density)   # id 接在骨架之后
            fines.append(fine)
            fine_bins.add(len(fines) - 1, x, y, z)
            break
        else:   # 6000 次都失败 → 孔隙里已经塞不下了
            raise RuntimeError(f"could only place {len(fines)} of {fine_count} fines without burying them in the skeleton")
    return fines


def evenly_spaced_centers(lower: float, upper: float, count: int) -> list[float]:
    """在 [lower, upper] 上等距取 count 个点（含两端）。"""
    if count < 2:
        raise ValueError("at least two screen wires are required")
    return [lower + i * (upper - lower) / (count - 1) for i in range(count)]


def build_screen(
    first_id: int,
    screen_x: float,
    ly: float,
    lz: float,
    wire_diameter: float,
    wire_count: int,
    density: float,
) -> tuple[list[tuple[int, int, float, float, float, float, float]], float, float]:
    """用固定的 type-3 珠子在 x=screen_x 平面上排出方孔编织筛网，返回 (珠子列表, 丝距, 名义孔径)。

    思路：筛网 = wire_count 根平行于 z 的丝 + wire_count 根平行于 y 的丝。
    每根丝用一串直径 = 丝径、彼此相切的珠子表示，这样 IBM 合并后沿丝方向没有缝隙。
    """
    radius = 0.5 * wire_diameter
    if wire_diameter <= 0.0 or wire_diameter >= min(ly, lz):
        raise ValueError("screen wire diameter is incompatible with the channel")
    # 丝的位置：在 y、z 方向各等距排 wire_count 根（贴边放，首尾丝的珠子与壁面相切）
    y_wires = evenly_spaced_centers(radius, ly - radius, wire_count)
    z_wires = evenly_spaced_centers(radius, lz - radius, wire_count)
    # 每根丝上的珠子个数：让相邻珠子中心距 <= 丝径（即相切或轻微重叠），避免轴向漏缝
    n_y = max(2, int(math.ceil((ly - wire_diameter) / wire_diameter)) + 1)
    n_z = max(2, int(math.ceil((lz - wire_diameter) / wire_diameter)) + 1)
    y_beads = evenly_spaced_centers(radius, ly - radius, n_y)
    z_beads = evenly_spaced_centers(radius, lz - radius, n_z)
    positions: set[tuple[float, float]] = set()   # 用集合去掉两组丝交叉点上的重复珠子
    for y in y_wires:          # 平行于 z 的丝：固定 y，沿 z 排珠子
        for z in z_beads:
            positions.add((y, z))
    for z in z_wires:          # 平行于 y 的丝：固定 z，沿 y 排珠子
        for y in y_beads:
            positions.add((y, z))
    screen = [
        (first_id + index, 3, radius, screen_x, y, z, density)
        for index, (y, z) in enumerate(sorted(positions))
    ]
    # 名义孔径 = 丝心距 − 丝径（相邻两根丝之间的净空）
    pitch = (ly - wire_diameter) / (wire_count - 1)
    aperture = pitch - wire_diameter
    if aperture <= 0.0:
        raise ValueError("screen wires overlap across the nominal aperture")
    return screen, pitch, aperture


def write_pack(path: Path, particles: list[tuple[int, int, float, float, float, float, float]], lx: float, ly: float, lz: float) -> None:
    """按 LIGGGHTS data 文件格式写出：头部(原子数、类型数、盒子)，Atoms 段每行 id type diameter density x y z。"""
    with path.open("w", encoding="ascii") as handle:
        handle.write("water_fangsha: Tang 2024 core-I target plus fixed outlet screen\n\n")
        handle.write(f"{len(particles)} atoms\n3 atom types\n\n")
        handle.write(f"0.0 {lx:.10e} xlo xhi\n0.0 {ly:.10e} ylo yhi\n0.0 {lz:.10e} zlo zhi\n\nAtoms\n\n")
        for pid, ptype, radius, x, y, z, density in particles:
            handle.write(f"{pid:d} {ptype:d} {2.0 * radius:.10e} {density:.6f} {x:.10e} {y:.10e} {z:.10e}\n")


def main() -> None:
    """解析参数 → 生成骨架、细颗粒、筛网 → 写文件并打印统计量。"""
    config = load_parameters()
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--seed", type=int, default=int(text(config, "PACK_RANDOM_SEED")))          # 固定种子保证可复现
    parser.add_argument("--lx-mm", type=float, default=real(config, "DOMAIN_LX_M") * 1000.0)
    parser.add_argument("--ly-mm", type=float, default=real(config, "DOMAIN_LY_M") * 1000.0)
    parser.add_argument("--lz-mm", type=float, default=real(config, "DOMAIN_LZ_M") * 1000.0)
    parser.add_argument("--core-x0-mm", type=float, default=real(config, "CORE_X0_M") * 1000.0)     # 岩心段起点
    parser.add_argument("--core-x1-mm", type=float, default=real(config, "CORE_X1_M") * 1000.0)     # 岩心段终点
    parser.add_argument("--target-porosity", type=float, default=real(config, "REFERENCE_CORE_POROSITY"))
    parser.add_argument("--fine-fraction", type=float, default=real(config, "FINE_FRACTION"))        # 细颗粒体积分数
    parser.add_argument("--max-skeleton-overlap", type=float, default=real(config, "SKELETON_MAX_OVERLAP_FRACTION"))
    parser.add_argument("--lattice-a-mm", type=float, default=real(config, "PACK_LATTICE_A_MM"))     # FCC 晶格常数
    parser.add_argument("--spacing-um", type=float, default=real(config, "REFERENCE_DIAMETER_UM") / real(config, "N"))  # 体素=LBM 格距
    parser.add_argument("--screen-x-mm", type=float, default=real(config, "SCREEN_X_M") * 1000.0)
    parser.add_argument("--screen-wire-diameter-um", type=float, default=real(config, "SCREEN_WIRE_DIAMETER_UM"))
    parser.add_argument("--screen-wire-count", type=int, default=int(text(config, "SCREEN_WIRE_COUNT")))
    parser.add_argument("--screen-density-kg-m3", type=float, default=real(config, "SCREEN_DENSITY_KG_M3"))
    parser.add_argument("--fine-diameter-um", type=float, default=real(config, "FINE_DIAMETER_UM"))
    parser.add_argument("--skeleton-density-kg-m3", type=float, default=real(config, "SKELETON_DENSITY_KG_M3"))
    parser.add_argument("--fine-density-kg-m3", type=float, default=real(config, "FINE_DENSITY_KG_M3"))
    parser.add_argument("--skeleton-diameter-min-mm", type=float, default=real(config, "SKELETON_DIAMETER_MIN_MM"))
    parser.add_argument("--skeleton-diameter-max-mm", type=float, default=real(config, "SKELETON_DIAMETER_MAX_MM"))
    parser.add_argument("--skeleton-diameter-mean-mm", type=float, default=real(config, "SKELETON_DIAMETER_MEAN_MM"))
    parser.add_argument("--skeleton-diameter-sigma-mm", type=float, default=real(config, "SKELETON_DIAMETER_SIGMA_MM"))
    args = parser.parse_args()
    # 命令行用 mm/µm 方便阅读，内部统一换成 m
    lx, ly, lz = [value * 1.0e-3 for value in (args.lx_mm, args.ly_mm, args.lz_mm)]
    core_x0, core_x1 = [value * 1.0e-3 for value in (args.core_x0_mm, args.core_x1_mm)]
    spacing = args.spacing_um * 1.0e-6
    screen_x = args.screen_x_mm * 1.0e-3
    wire_diameter = args.screen_wire_diameter_um * 1.0e-6
    lattice_a = args.lattice_a_mm * 1e-3
    if not (0.0 < core_x0 < core_x1 <= lx and 0.0 < args.fine_fraction < 1.0 and 0.0 < args.target_porosity < 1.0 and 0.0 <= args.max_skeleton_overlap < 0.5 and lattice_a > 0.0):
        raise ValueError("invalid geometry or packing parameters")
    if not (core_x1 < screen_x < lx and args.screen_wire_count >= 2 and args.screen_density_kg_m3 > 0.0):
        raise ValueError("invalid outlet screen parameters")   # 筛网必须在岩心下游、出口上游
    rng = random.Random(args.seed)
    # ── 三步生成：骨架 → 细颗粒 → 筛网；id 依次连续编号 ──
    skeleton, solid_mask = build_skeleton(
        rng, core_x0, core_x1, ly, lz, spacing, args.target_porosity, args.max_skeleton_overlap, lattice_a,
        args.skeleton_diameter_min_mm * 1e-3, args.skeleton_diameter_max_mm * 1e-3,
        args.skeleton_diameter_mean_mm * 1e-3, args.skeleton_diameter_sigma_mm * 1e-3,
        args.skeleton_density_kg_m3,
    )
    fines = build_fines(
        rng, skeleton, solid_mask, spacing, core_x0, lx, ly, lz, args.fine_fraction,
        args.fine_diameter_um * 1e-6, args.fine_density_kg_m3, args.skeleton_diameter_max_mm * 1e-3,
    )
    screen, pitch, aperture = build_screen(
        len(skeleton) + len(fines) + 1,
        screen_x,
        ly,
        lz,
        wire_diameter,
        args.screen_wire_count,
        args.screen_density_kg_m3,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_pack(args.output, skeleton + fines + screen, lx, ly, lz)
    # ── 打印统计量，便于和 validate_water_fangsha_pack.py 的结果对照 ──
    solid_fraction = float(np.mean(solid_mask))
    fine_volume = len(fines) * 4.0 * math.pi * (0.5 * args.fine_diameter_um * 1e-6) ** 3 / 3.0
    skeleton_volume = float(np.count_nonzero(solid_mask)) * spacing**3
    print(f"skeleton_count={len(skeleton)}")
    print(f"skeleton_union_porosity={1.0 - solid_fraction:.8f}")
    print(f"fine_count={len(fines)}")
    print(f"fine_concentration_vs_skeleton_union={fine_volume / (fine_volume + skeleton_volume):.8f}")
    print(f"screen_count={len(screen)}")
    print(f"screen_x_mm={screen_x * 1e3:.8f}")
    print(f"screen_pitch_um={pitch * 1e6:.8f}")
    print(f"screen_nominal_aperture_um={aperture * 1e6:.8f}")
    print(f"screen_aperture_to_fine_ratio={aperture / (args.fine_diameter_um * 1e-6):.8f}")   # >1 细颗粒才可能穿网
    print(f"output={args.output}")


if __name__ == "__main__":
    main()

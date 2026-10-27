#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_liggghts_to_csv.py —— 把 LIGGGHTS 颗粒 dump 文件(particles_*.liggghts)转成 CSV。

【这个脚本做什么】
  LIGGGHTS 的 dump 是自定义文本格式，Excel/pandas/ParaView 不能直接读。
  本脚本把一个 dump 里的所有帧展开成“每颗粒一行”的 CSV，并把时间步和盒子边界
  作为列附在每一行前面，方便后续筛选和绘图。

【什么时候、由谁调用】
  - 自动：C++ 主程序 oilwaterDrive 结束时，对 post/ 下每个 <前缀>particles_*.liggghts
    各调用一次本脚本；任一次失败主程序以退出码 5 结束。
  - 手工：python3 convert_liggghts_to_csv.py <某个 particles_xxx.liggghts>

【输入格式】（每一“帧”= 某个时间步的全部颗粒快照，一个文件可含多帧）
    ITEM: TIMESTEP
    12345                      <- 该帧对应的 DEM 步数
    ITEM: NUMBER OF ATOMS
    502                        <- 该帧颗粒个数
    ITEM: BOX BOUNDS pp pp pp
    0.0 2.3e-3                 <- x 方向盒子下界 上界
    0.0 1.0e-3                 <- y 方向
    0.0 1.0e-3                 <- z 方向
    ITEM: ATOMS id type x y z vx vy vz radius fx fy fz   <- 列名
    1 1 ... (502 行颗粒数据)

【输出】
  与输入同目录、同名、后缀改为 .csv。
  列：timestep, box_xlo, box_xhi, box_ylo, box_yhi, box_zlo, box_zhi, 然后是 dump 的原始各列。

【判据】
  头部标记缺失或文件被截断（颗粒行数不够）时抛 ValueError，不写出半截 CSV。
"""
from __future__ import annotations   # 让类型注解(如 -> Path)在旧版 Python 也能用
import argparse                       # 解析命令行参数
import csv                            # 标准库 CSV 写出器(自动处理分隔符/换行)
from pathlib import Path             # 面向对象的路径操作(比字符串拼接更安全)

def parse_one_dump(path: Path):
    """读入一个 dump 文件, 解析出其中所有帧, 返回帧列表。"""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()  # 整文件读入并按行切分
    i = 0                             # 当前扫描到的行号(从 0 开始)
    frames = []                       # 收集所有帧的容器

    while i < len(lines):             # 从头到尾扫描每一行
        if not lines[i].startswith("ITEM: TIMESTEP"):  # 没遇到帧起始标记
            i += 1                    # 跳过该行, 继续找下一帧的开头
            continue                  # 回到 while 顶部

        timestep = int(lines[i+1].strip())  # 帧标记的下一行就是 DEM 步数

        if not lines[i+2].startswith("ITEM: NUMBER OF ATOMS"):  # 校验文件格式(防止读错)
            raise ValueError(f"{path}: NUMBER OF ATOMS header missing near line {i+3}")
        n_atoms = int(lines[i+3].strip())   # 该帧颗粒个数

        if not lines[i+4].startswith("ITEM: BOX BOUNDS"):  # 校验盒子边界段存在
            raise ValueError(f"{path}: BOX BOUNDS header missing near line {i+5}")

        xlo, xhi = map(float, lines[i+5].split()[:2])  # x 盒子下/上界(取前两个数)
        ylo, yhi = map(float, lines[i+6].split()[:2])  # y 盒子下/上界
        zlo, zhi = map(float, lines[i+7].split()[:2])  # z 盒子下/上界

        atom_header = lines[i+8].strip()    # 第 8 行是列名行 "ITEM: ATOMS id type x ..."
        if not atom_header.startswith("ITEM: ATOMS"):  # 校验列名行存在
            raise ValueError(f"{path}: ATOMS header missing near line {i+9}")

        cols = atom_header.replace("ITEM: ATOMS", "", 1).strip().split()  # 去掉前缀, 拆出每个列名
        start = i + 9                  # 颗粒数据的第一行行号
        end = start + n_atoms          # 颗粒数据的结束行号(不含)
        if end > len(lines):           # 文件被截断(写一半中断)
            raise ValueError(f"{path}: unexpected EOF while reading atoms")

        frames.append({                # 把这一帧的信息打包存起来
            "timestep": timestep,      # 该帧步数
            "box": (xlo, xhi, ylo, yhi, zlo, zhi),  # 盒子边界(6 个数)
            "cols": cols,              # 列名列表
            "rows": [lines[k].split() for k in range(start, end)]  # 每个颗粒一行(已拆成字段列表)
        })
        i = end                        # 跳到本帧末尾, 继续找下一帧

    return frames                      # 返回解析出的所有帧

def convert_file(path: Path):
    """把单个 dump 文件转换成同名 .csv, 返回输出路径。"""
    frames = parse_one_dump(path)      # 先解析出所有帧
    if not frames:                     # 一帧都没有(空文件/格式错)
        raise ValueError(f"{path}: no dump frames found")

    out = path.with_suffix(".csv")     # 输出文件名: 把后缀换成 .csv
    with out.open("w", newline="", encoding="utf-8") as f:  # newline="" 是 csv 模块推荐写法
        writer = csv.writer(f)         # 创建 CSV 写出器
        header = [                     # CSV 表头: 时间步 + 盒子边界 + 颗粒各列
            "timestep",
            "box_xlo", "box_xhi",
            "box_ylo", "box_yhi",
            "box_zlo", "box_zhi",
            *frames[0]["cols"]         # 展开第一帧的列名(各帧列名一致)
        ]
        writer.writerow(header)        # 写表头(只写一次)

        for fr in frames:              # 逐帧写出
            xlo, xhi, ylo, yhi, zlo, zhi = fr["box"]  # 取该帧盒子边界
            for row in fr["rows"]:     # 该帧每个颗粒一行
                writer.writerow([      # 每行 = 时间步 + 盒子边界 + 该颗粒字段
                    fr["timestep"],
                    xlo, xhi, ylo, yhi, zlo, zhi,
                    *row               # 展开该颗粒的所有字段
                ])
    return out                         # 返回 CSV 路径

def main():
    """命令行入口: 读取一个文件路径参数, 转换并打印结果。"""
    parser = argparse.ArgumentParser()  # 创建参数解析器
    parser.add_argument("path")         # 唯一位置参数: dump 文件路径
    args = parser.parse_args()          # 解析命令行
    path = Path(args.path).resolve()    # 转成绝对路径(避免相对路径歧义)
    out = convert_file(path)            # 执行转换
    print(f"[OK] {path.name} -> {out.name}")  # 打印 "输入名 -> 输出名"

if __name__ == "__main__":             # 只有被直接运行时才执行(被 import 时不执行)
    main()

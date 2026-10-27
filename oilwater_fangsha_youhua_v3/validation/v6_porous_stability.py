# -*- coding: utf-8 -*-
"""
验证 6：二维多孔介质中的稳定性（2026-10-26 新增）
------------------------------------------------------------------------
为什么要有这个验证：
  v1–v5 都是宽通道/平板，喉道很宽。真实算例（N=20 时）颗粒间的喉道只有 1–4 格宽，
  壁面润湿力在窄喉道里两侧同时作用，最容易把流体挤爆。本机测试第 7 步发散就是这样被发现的。
设置：200×90 的二维区域里随机放 27 个直径 20–37 格的"骨架"圆盘（部分重叠，孔隙度约 0.24）
      和 25 个贴在骨架表面的 5 格"细粒"圆盘；整个孔隙先充满油，左侧 12 格是水，
      加一个很小的体积力往右推。和 C++ 一样用压力平衡的初始界面（set_composition_pb）。
输出：每隔一段打印孔隙里总密度的最小/最大值、最大混合物速度；出现 NaN 即判为发散。
2026-10-26 结果（τ_w=0.7、τ_o=5.45 即 M≈20，G=2.5，几何接触角，WETTING_GAIN=2.4，3000 步）：
      θ = 30°、45°、90°、135°、150° 全部稳定，孔隙内总密度始终在 0.57–1.15 之间。
      作为对比，旧的虚拟密度写法（LBM2P_WALL=repel，s_wall=0.3 或 0.7）在 500 步内发散。
2026-10-27 加入修订八后复测：θ=45° 总密度 [0.570, 1.078]，θ=135° [0.572, 1.044]，3000 步稳定。
用法：python3 v6_porous_stability.py [θ（度）] [步数]
      例：python3 v6_porous_stability.py 45 3000
"""
import os
import sys
import numpy as np
from lbm2p import TwoPhase

theta = float(sys.argv[1]) if len(sys.argv) > 1 else 90.0
nsteps = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
TW, TO = float(os.environ.get('TW', '0.7')), float(os.environ.get('TO', '5.4505'))

rng = np.random.default_rng(1)
NX, NY = 200, 90
x = np.arange(NX)[:, None]; y = np.arange(NY)[None, :]
solid = np.zeros((NX, NY), bool)
cents = []
for k in range(400):
    r = rng.uniform(10, 18.5); cx = rng.uniform(25, NX-5); cy = rng.uniform(0, NY)
    if any((cx-a)**2 + (cy-b)**2 < (0.8*(r+c))**2 for a, b, c in cents):
        continue
    cents.append((cx, cy, r))
for cx, cy, r in cents:
    for dy in (-NY, 0, NY):
        solid |= (x-cx)**2 + (y-cy-dy)**2 < r*r
for cx, cy, r in cents[:25]:      # 细粒：贴在骨架表面
    ang = rng.uniform(0, 2*np.pi); fx = cx + (r+2.0)*np.cos(ang); fy = (cy + (r+2.0)*np.sin(ang)) % NY
    solid |= (x-fx)**2 + (y-fy)**2 < 2.5**2
porosity = 1 - solid[25:].mean()

sim = TwoPhase((NX, NY, 1), tau_w=TW, tau_o=TO, G=2.5, theta=theta, body_force=(2e-6, 0, 0))
sim.set_solid(solid[:, :, None])
sim.set_composition_pb((x < 12)[:, :, None] & ~solid[:, :, None], width=2)
print(f"θ={theta}° τ={TW}/{TO} 孔隙度={porosity:.3f} 骨架圆盘数={len(cents)} 壁面模型={sim.wall_model}", flush=True)
fl = ~solid
for t in range(1, nsteps + 1):
    sim.step()
    if t in (1, 10, 50) or t % 500 == 0:
        rw, ro, u, p = sim.fields()
        rt = (rw + ro)[:, :, 0][fl]
        um = np.sqrt((u[:, :, :, 0]**2).sum(0))[fl]
        bad = not np.isfinite(rt).all()
        print(f"  第 {t:5d} 步：孔隙内总密度 [{np.nanmin(rt):.3f}, {np.nanmax(rt):.3f}]  最大速度 {np.nanmax(um):.3e}"
              f"{'  ← NaN，发散' if bad else ''}", flush=True)
        if bad:
            break

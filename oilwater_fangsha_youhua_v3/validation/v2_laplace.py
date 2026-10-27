# -*- coding: utf-8 -*-
"""
验证 2：Laplace 定律 —— 标定界面张力 σ
--------------------------------------
原理：静止液滴内外压差 Δp = σ/R（二维圆滴）。不同半径各算一次，Δp 对 1/R 作图
     应是过原点的直线，斜率就是格子单位的界面张力 σ_lb。
换算：σ_phys = σ_lb · ρ_ref · Δx³ / Δt²。Δt 由 τ_w 决定：Δt = (τ_w−0.5)/3 · Δx² / ν_w。
     因此 σ、τ、Δx 三者绑在一起，不能各自随便选——这就是为什么方法节要报告
     毛细数 Ca = μ u / σ 而不是只报告 σ。
用法：python3 v2_laplace.py [G]   （默认 G=2.5，τ_w=0.8，τ_o=2.5）
"""
import sys
import numpy as np
from lbm2p import TwoPhase

G = float(sys.argv[1]) if len(sys.argv) > 1 else 2.5
import os
# 2026-10-25：τ 可用环境变量改，例如  TW=0.6 TO=2.975 python3 v2_laplace.py
TAU_W, TAU_O = float(os.environ.get('TW', '0.8')), float(os.environ.get('TO', '2.5'))
N = 64


def droplet(R, nsteps=5000):
    sim = TwoPhase((N, N, 1), tau_w=TAU_W, tau_o=TAU_O, G=G)
    x = np.arange(N)[:, None, None]; y = np.arange(N)[None, :, None]
    sim.set_composition((x - N/2)**2 + (y - N/2)**2 < R**2)   # 中央水滴，周围油
    ok = sim.run(nsteps)
    rw, ro, u, p = sim.fields()
    c = N//2
    dp = p[c-2:c+3, c-2:c+3, 0].mean() - p[:4, :4, 0].mean()
    Reff = np.sqrt((rw[:, :, 0] > 0.5*(rw.max() + rw.min())).sum()/np.pi)
    return ok, Reff, dp, np.abs(u).max()


rows = []
for R in (8, 11, 14, 18):
    ok, Reff, dp, umax = droplet(R)
    print(f"R={R:2d}: {'ok' if ok else '发散'} R_eff={Reff:.2f} dp={dp:.5f} dp*R={dp*Reff:.4f} 伪速度={umax:.1e}", flush=True)
    if ok:
        rows.append((Reff, dp))
X = np.array([1/r for r, _ in rows]); Y = np.array([d for _, d in rows])
s, c0 = np.linalg.lstsq(np.vstack([X, np.ones_like(X)]).T, Y, rcond=None)[0]
print(f"\nσ_lb = {s:.4f}（截距 {c0:.1e}，应接近 0）")
# 用生产参数做物理换算：Δx=5 μm，ν_w=1.01e-6 m²/s，ρ_ref=1000 kg/m³
dx, nu_w, rho_ref = 5e-6, 1.01e-6, 1000.0
dt = (TAU_W - 0.5)/3*dx**2/nu_w
print(f"生产参数下 Δt = {dt:.3e} s，σ_phys = {s*rho_ref*dx**3/dt**2*1e3:.2f} mN/m（真实油水约 20–30 mN/m）")
print("→ 若 σ_phys 偏小，应通过降低入口速度来匹配毛细数 Ca，而不是把 G 调到失稳。")

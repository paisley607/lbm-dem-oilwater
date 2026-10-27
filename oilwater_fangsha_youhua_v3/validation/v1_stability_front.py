# -*- coding: utf-8 -*-
"""
验证 1：一维油水前缘的稳定性（对应审计报告附录 A 算例 1 的"修复后"版本）
----------------------------------------------------------------------------
目的：确认新算法（公共速度 + 精确力注入 + 局部混合松弛时间）在给定 τ_w、τ_o、G 下
     能稳定分相，并给出两相共存密度（主相/次相）。
2026-10-25 结果（G=2.5，混合 τ）：0.6/2.56 稳定（水区溶油 0.15%）；0.8/7.925 稳定（0.5%）。
     旧做法（LBM2P_FIXED_TAU=1）在 0.6/3.0 时水区溶油 5%，水的有效粘度被"稀释"成 2 倍以上。
判据：跑 20000 步不出现 NaN；水区 ρ_w≈1、ρ_o≈残余；油区反之。
用法：python3 v1_stability_front.py [tau_w] [tau_o] [G]
"""
import sys
import numpy as np
from lbm2p import TwoPhase

tau_w = float(sys.argv[1]) if len(sys.argv) > 1 else 0.8
tau_o = float(sys.argv[2]) if len(sys.argv) > 2 else 2.5
G = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0

nx = 120
sim = TwoPhase((nx, 1, 1), tau_w=tau_w, tau_o=tau_o, G=G)
x = np.arange(nx)
sim.set_composition(((x >= nx//4) & (x < 3*nx//4))[:, None, None])   # 1 格阶跃，故意用最苛刻的初值
ok = sim.run(20000)
rw, ro, u, p = sim.fields()
print(f"tau_w={tau_w} tau_o={tau_o} G={G}: {'稳定' if ok else '发散'}")
if ok:
    print(f"  水区: rho_w={rw[nx//2,0,0]:.4f} rho_o={ro[nx//2,0,0]:.4f}")
    print(f"  油区: rho_w={rw[0,0,0]:.4f} rho_o={ro[0,0,0]:.4f}")
    print(f"  最大伪速度 |u| = {np.abs(u).max():.2e}（越小越好，典型 1e-3 量级）")

# -*- coding: utf-8 -*-
"""
验证 3：平板接触角 —— 检查/标定润湿参数 WALL_CONTACT_ANGLE_DEG 与 WETTING_GAIN
------------------------------------------------------------------------
原理（2026-10-26 起与 C++ 一致，用"几何接触角"）：
      在平壁上放一个半圆水滴，周围是油。给壁面一个"输入接触角" θ_in，
      程序在接触线附近按 θ_in 外推固体邻格的相场（见 lbm2p.py 说明），稳定后
      用球冠几何从液滴高度 h 和底宽 2a 反算"实测接触角"：R=(a²+h²)/(2h)，cosθ=(R−h)/R。
      输入角和实测角越接近，说明 WETTING_GAIN 标定得越好。
2026-10-26 结果（τ_w=0.7、τ_o=5.45 即 M≈20，G=2.5，8000 步）：
      WETTING_GAIN=1.0： 30→51°   60→76°   90→88°   120→102°   150→125°（不够"强"）
      WETTING_GAIN=2.4： 30→25°   60→58°   90→88°   120→118°            ← 默认值
2026-10-27 加入修订八（痕量组分松弛过渡 + 密度地板，见 lbm2p.py）后：45→40.5°  90→88.4°  135→129.0°
      （油湿一侧偏小约 6°；需要精确角度时可改 GAIN 重新标定）
旧写法（虚拟密度 s_wall，修订三）在二维多孔介质里 s=0.3、0.7 都发散，已不再使用；
可用 LBM2P_WALL=repel 切回旧写法对比（此时参数是 s_wall，不是角度）。
用法：python3 v3_contact_angle.py                 （扫 30 60 90 120 150 度）
      python3 v3_contact_angle.py 45 135          （只算指定角度）
      GAIN=2.0 TW=0.8 TO=2.5 python3 v3_contact_angle.py 60   （换系数、换 τ 重新标定）
"""
import os
import sys
import numpy as np
from lbm2p import TwoPhase

NX, NY, WALL = 120, 60, 4
G = 2.5
TAU_W, TAU_O = float(os.environ.get('TW', '0.7')), float(os.environ.get('TO', '5.4505'))
GAIN = float(os.environ.get('GAIN', '2.4'))
NSTEPS = int(os.environ.get('NSTEPS', '8000'))


def contact_angle(theta_in, R0=18, nsteps=NSTEPS):
    sim = TwoPhase((NX, NY, 1), tau_w=TAU_W, tau_o=TAU_O, G=G, theta=theta_in, wet_gain=GAIN)
    x = np.arange(NX)[:, None, None]; y = np.arange(NY)[None, :, None]
    solid = np.zeros((NX, NY, 1), bool); solid[:, :WALL] = True
    sim.set_solid(solid)
    sim.set_composition(((x - NX/2)**2 + (y - WALL)**2 < R0**2) & ~solid)
    ok = sim.run(nsteps)
    rw, ro, u, p = sim.fields()
    phi = ((rw - ro)/(rw + ro))[:, :, 0]
    col = (phi > 0) & ~solid[:, :, 0]
    if not ok or not col.any():
        return None
    ys = np.where(col.any(axis=0))[0]
    h = ys.max() - WALL + 1
    base = col[:, WALL].sum()
    if base < 2:
        return 180.0
    a = base/2.0
    Rc = (a*a + h*h)/(2*h)
    return float(np.degrees(np.arccos(np.clip((Rc - h)/Rc, -1, 1))))


if __name__ == "__main__":
    vals = [float(v) for v in sys.argv[1:]] or [30, 60, 90, 120, 150]
    print(f"τ_w={TAU_W} τ_o={TAU_O} G={G} WETTING_GAIN={GAIN}")
    print("输入θ(度)  实测θ(度)  含义")
    for t in vals:
        th = contact_angle(t)
        desc = "—" if th is None else ("水湿" if th < 75 else ("中性" if th <= 105 else "油湿"))
        print(f"{t:7.1f}    {('发散' if th is None else f'{th:7.1f}')}    {desc}", flush=True)

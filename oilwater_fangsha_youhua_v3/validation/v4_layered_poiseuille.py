# -*- coding: utf-8 -*-
"""
验证 4：两相分层 Poiseuille 流 —— 检查粘度比和两相之间的粘性耦合
------------------------------------------------------------------
设置：两块平壁之间，下半是水、上半是油，x 方向周期并加体积力 g 驱动。
解析解：两层流体各自是抛物线，在界面处速度连续、剪应力连续：
      μ_w u_w'(0) = μ_o u_o'(0)。
输出（2026-10-25 改为更有意义的指标）：
  · L2 相对误差：整条剖面与解析解的均方根偏差（旧版用"单点最大误差"，会被界面上 1 个格点的
    尖峰主导，看不出整体好坏）；
  · 总流量误差、水层/油层各自的流量误差：渗透率、相对渗透率都是由流量算的，这个最关键。
2026-10-25 结果（G=2.5，混合 τ + 只限偏离量，详见 修改说明_20261025.md）：
     τ_w/τ_o      运动粘度比   L2 误差   总流量   水层流量   油层流量
     0.8/0.8         1         4.5%     +1.7%     +0.8%     +2.6%
     0.8/2.5         6.7       4.4%     +4.0%     +0.2%    +11.2%
     0.8/7.925      24.8      10.1%     +3.4%     −3.7%    +35.6%
     0.6/1.12        6.2       5.3%     +2.2%     −1.5%     +9.1%
     0.6/2.56       20.6       9.9%     −1.1%     −6.0%    +17.8%
     0.7/0.995       2.5       6.0%     −2.5%     −3.8%     −0.8%
     0.7/5.45       24.8       8.0%     +1.0%     −3.4%    +20.9%
     0.8/19.0       61.7      13.0%     +1.8%     −5.4%    +72.5%
  油层（粘的那一层）流量偏大，来自 Shan-Chen 在壁面旁的"密度亏损层"造成的表观滑移，
  粘度比越大越明显——这是该模型的已知局限，方法节里应当说明。
  旧做法（LBM2P_FIXED_TAU=1，整体限速）：0.8/0.8 时界面变成"墙"，L2 误差约 90%；0.6/1.12 发散。
★★ 2026-10-26：组分间动量交换改用 τ_D=1（DIFFUSION_TAU），壁面改为几何接触角（θ=90° 时壁面不施力）。
  τ 0.7/5.45 的结果见 README 验证表（10-25 版油层流量 +20.9%，τ_D=1 后约 +3%）。
  固定 τ（mix=False）算例仍按旧做法（τ_D=0）跑，用于 v5 的新旧对比。
★★ 2026-10-27：加入修订八（痕量组分松弛过渡 + 密度地板）后，τ 0.7/5.45：L2 2.1%，总流量 +0.1%，水层 −0.6%，油层 +3.2%
  （同一版代码 LBM2P_TRACE=0 关掉修订八：L2 1.7%，+0.7%，+0.3%，+2.7%）。
用法：python3 v4_layered_poiseuille.py [τ_w τ_o] [fixed]
      例：python3 v4_layered_poiseuille.py 0.6 2.975        （默认配置 M≈20 对应的 τ）
"""
import sys
import numpy as np
import numpy.linalg as la
from lbm2p import TwoPhase

NY, WALL = 64, 2            # 流体区 y = 2 … 61（60 格），界面在中间
G = 2.5
GX = 1e-6


def analytic(nu_w, nu_o, H, yy):
    """两层 Poiseuille 解析解。未知数 [Aw, Bw, Ao, Bo]：壁面无滑移×2、界面速度连续、剪应力连续。"""
    M = np.array([[-H, 1, 0, 0],           # u_w(−H)=0
                  [0, 0, H, 1],            # u_o(+H)=0
                  [0, 1, 0, -1],           # u_w(0)=u_o(0)
                  [nu_w, 0, -nu_o, 0]])    # μ_w u_w'(0) = μ_o u_o'(0)（两层格子密度都≈1）
    rhs = np.array([GX/(2*nu_w)*H*H, GX/(2*nu_o)*H*H, 0, 0])
    Aw, Bw, Ao, Bo = la.solve(M, rhs)
    return np.where(yy < 0, -GX/(2*nu_w)*yy**2 + Aw*yy + Bw, -GX/(2*nu_o)*yy**2 + Ao*yy + Bo)


def run_case(tau_w, tau_o, mix=True, nsteps=None, verbose=True):
    """跑一个分层流算例，返回指标字典（发散时 stable=False）。"""
    if nsteps is None:                       # τ_w 越小粘度越小、收敛越慢
        nsteps = 40000 if tau_w >= 0.7 else 60000
    sim = TwoPhase((1, NY, 1), tau_w=tau_w, tau_o=tau_o, G=G, theta=90.0, body_force=(GX, 0, 0), mix_tau=mix,
                   tau_d=(1.0 if mix else 0.0))      # 固定 τ 对照组保持旧做法
    y = np.arange(NY)[None, :, None]
    solid = np.zeros((1, NY, 1), bool); solid[:, :WALL] = True; solid[:, NY-WALL:] = True
    sim.set_solid(solid)
    sim.set_composition((y < NY//2) & ~solid)
    ok = sim.run(nsteps)
    out = dict(tau_w=tau_w, tau_o=tau_o, mix=sim.mix_tau, stable=ok)
    label = '混合τ' if sim.mix_tau else '固定τ'
    if not ok:
        if verbose:
            print(f"τ_w={tau_w} τ_o={tau_o} {label}：发散")
        return out
    rw, ro, u, _ = sim.fields()
    ux, rw0, ro0 = u[0, 0, :, 0], rw[0, :, 0], ro[0, :, 0]
    H = (NY - 2*WALL)/2.0
    yy = np.arange(NY) - (WALL + H - 0.5)            # 以界面为原点
    nu_w, nu_o = (tau_w - 0.5)/3, (tau_o - 0.5)/3
    ua = analytic(nu_w, nu_o, H, yy)
    fl = ~solid[0, :, 0]
    lower = fl & (np.arange(NY) < NY//2)
    upper = fl & (np.arange(NY) >= NY//2)
    out.update(
        L2=float(np.sqrt(((ux - ua)[fl]**2).sum()/(ua[fl]**2).sum())),
        q_total=float(((rw0 + ro0)*ux)[fl].sum()/ua[fl].sum() - 1),
        q_water=float((rw0*ux)[fl].sum()/ua[lower].sum() - 1),
        q_oil=float((ro0*ux)[fl].sum()/ua[upper].sum() - 1),
        water_layer=(float(rw0[WALL + int(H/2)]), float(ro0[WALL + int(H/2)])),
        oil_layer=(float(rw0[NY - WALL - int(H/2)]), float(ro0[NY - WALL - int(H/2)])),
        profile=(ux, ua))
    if verbose:
        print(f"τ_w={tau_w} τ_o={tau_o}（运动粘度比 {nu_o/nu_w:.1f}）{label}：稳定")
        print(f"  L2 相对误差 {out['L2']*100:.1f}%   总流量误差 {out['q_total']*100:+.1f}%   "
              f"水层流量 {out['q_water']*100:+.1f}%   油层流量 {out['q_oil']*100:+.1f}%")
        print(f"  水层组成 ρ_w={out['water_layer'][0]:.3f} ρ_o={out['water_layer'][1]:.3f}；"
              f"油层组成 ρ_w={out['oil_layer'][0]:.3f} ρ_o={out['oil_layer'][1]:.3f}")
    return out


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "fixed"]
    tw = float(args[0]) if len(args) > 0 else 0.8
    to = float(args[1]) if len(args) > 1 else 2.5
    r = run_case(tw, to, mix=("fixed" not in sys.argv))
    if r["stable"]:
        ux, ua = r["profile"]
        print("  y    模拟 u      解析")
        for k in range(WALL, NY - WALL, 6):
            print(f"{k:3d}  {ux[k]:.3e}  {ua[k]:.3e}")

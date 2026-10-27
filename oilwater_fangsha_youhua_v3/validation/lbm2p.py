# -*- coding: utf-8 -*-
"""
lbm2p.py —— 油水两相 Shan-Chen 格子 Boltzmann 的 numpy 小程序（验证算例公用库）
=====================================================================================

这个文件的作用：
  用纯 Python（numpy）把 C++ 主程序 oilwaterShanChen3D.hh 里的两相算法原样复现一遍，
  拿来跑小型验证算例（Laplace 定律、接触角、分层 Poiseuille 流、稳定性）。
  它和 C++ 用的是同一套公式，所以在这里验证通过的参数（G、τ、壁面润湿参数），
  在 C++ 里也应该是合理的。它不是生产求解器，网格很小、没有并行、没有颗粒运动。

先补几个基础概念（初学者请先读这一段）：
  * 格子 Boltzmann 方法（LBM）：不直接解 Navier-Stokes 方程，而是在每个格点上
    存 19 个方向的"粒子分布函数" f_i（D3Q19），每一步做两件事：
      碰撞（collide）——让 f_i 向平衡态 feq 靠拢；
      迁移（stream）——把 f_i 沿各自方向搬到相邻格点。
    密度 ρ = Σ f_i，动量 j = Σ f_i c_i，速度 u = j/ρ。
  * τ（松弛时间）：决定流体粘度，ν_lb = (τ − 0.5)/3。τ 越接近 0.5 粘度越小、也越不稳定。
    多相模型对 τ 很敏感。2026-10-25 的验证：用"局部混合松弛时间"（mix_tau=True）时 τ_w=0.7 在分层流/液滴/接触角
    算例都稳定；τ_w=0.6 时小液滴和接触角算例发散；0.55 连分层流都发散。
  * Shan-Chen 多组分模型：油、水各一套 f_i。两种组分之间加一个"排斥力"
      F_w = −G ρ_w Σ_i t_i ρ_o(x+c_i) c_i   （水受油的推挤，推向油少的方向）
      F_o = −G ρ_o Σ_i t_i ρ_w(x+c_i) c_i
    G 越大排斥越强、界面张力越大；G 太大会失稳。
  * 公共速度 u'：两种组分碰撞时都向同一个速度 u' = Σω_σ j_σ / Σω_σ ρ_σ 靠拢，
    这是两相之间粘性耦合（切向动量交换）的来源。旧代码缺了这一步，所以会发散。
  * 润湿性（接触角）：求和 Σ_i t_i ρ(x+c_i) c_i 碰到固体邻格时，要给固体一个"虚拟密度"。
    ★★ 2026-10-26 起与 C++ 一致，默认用"几何接触角"（wall_model='geom'）：
      相场 φ=(ρ_w−ρ_o)/(ρ_w+ρ_o)，固体邻格的 φ 按接触角 θ 从本格外推
        φ_s = φ + c_i·∇_tφ − K|∇_tφ|(c_i·n)，K = gain·cotθ（gain 标定为 2.4），
      总密度取本格的值（镜像）。本体流体里 ∇φ≈0，壁面不施任何力，只有接触线附近起作用。
    旧写法（'repel'：只排斥非润湿相；'virtual'：ρ̃_w=s·ρ_wall、ρ̃_o=(1−s)·ρ_wall）保留做对比，
    它们会在整个壁面推挤非润湿相，窄喉道里会发散。
  * 组分间动量交换（2026-10-26）：u' 和 j* 里的 ω 用单独的 ω_D=1/τ_D（默认 τ_D=1），
    粘性仍按混合 τ 在碰撞里松弛。旧做法（τ_D=0，用粘性 ω）在 τ_o 大时界面动量不守恒。
  * 痕量组分（2026-10-27，C++ 修订八）：质量分数 < s_lo 的组分，粘性松弛频率从混合 ω 过渡到 ω_D；
    痕量组分密度 < 0.5·ρ_res 时重置为平衡态。本体组分完全不受影响。

本文件提供：
  TwoPhase 类：建立网格、初始化油水分布、设置固体格点、推进时间步、提取结果。
  pressure()：Shan-Chen 模型的总压力 p = c_s²(ρ_w+ρ_o) + c_s² G ρ_w ρ_o。
"""
import numpy as np

# ---------- D3Q19 格子常数 ----------
# 19 个离散速度方向 c_i（第 0 个是静止方向）
C = np.array([[0, 0, 0],
              [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
              [1, 1, 0], [-1, -1, 0], [1, -1, 0], [-1, 1, 0],
              [1, 0, 1], [-1, 0, -1], [1, 0, -1], [-1, 0, 1],
              [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, 1]], dtype=int)
W = np.array([1/3] + [1/18]*6 + [1/36]*12)        # 每个方向的权重 t_i
Q = 19                                             # 方向数
CS2 = 1.0/3.0                                      # 格子声速平方 c_s²
Cf = C.astype(float)
OPP = np.array([int(np.where((C == -C[i]).all(axis=1))[0][0]) for i in range(Q)])  # 相反方向的编号
# 用于"正则化碰撞"的二阶 Hermite 张量 Q_i = c_i c_i − c_s² I（存 6 个独立分量）
PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
MULT = np.array([1, 1, 1, 2, 2, 2], float)         # 对称张量非对角分量出现两次
QT = np.array([[Cf[i, a]*Cf[i, b] - (CS2 if a == b else 0.0) for (a, b) in PAIRS] for i in range(Q)])
PAIR_CC = np.array([[Cf[i, a]*Cf[i, b] for (a, b) in PAIRS] for i in range(Q)])


def feq(rho, u):
    """二阶平衡态分布：f_eq_i = t_i ρ [1 + 3 c·u + 4.5 (c·u)² − 1.5 u²]。"""
    cu = np.einsum('qd,d...->q...', Cf, u)
    usq = (u**2).sum(axis=0)
    return W[:, None, None, None]*rho[None]*(1 + 3*cu + 4.5*cu**2 - 1.5*usq[None])


def moments(f):
    """由分布函数求密度 ρ 和动量 j。"""
    rho = f.sum(axis=0)
    j = np.einsum('qd,q...->d...', Cf, f)
    return rho, j


def stream(f):
    """迁移步：把 f_i 沿方向 c_i 搬一格（周期边界，固体格点另行处理）。"""
    out = np.empty_like(f)
    for i in range(Q):
        out[i] = np.roll(f[i], shift=tuple(C[i]), axis=(0, 1, 2))
    return out


def neighbour_sum(rho):
    """H(x) = Σ_i t_i ρ(x+c_i) c_i —— 密度梯度型矢量，Shan-Chen 力的核心。"""
    H = np.zeros((3,) + rho.shape)
    for i in range(1, Q):
        nb = np.roll(rho, shift=tuple(-C[i]), axis=(0, 1, 2))
        for d in range(3):
            if C[i, d] != 0:
                H[d] += W[i]*nb*C[i, d]
    return H


def collide_regularized(f, omega):
    """正则化 BGK 碰撞（与 Palabos 的 RegularizedBGK 思路一致，这里不加 Smagorinsky）。
    只保留非平衡部分的二阶矩 Π_neq，再按 (1−ω) 衰减；比普通 BGK 更稳定。"""
    rho, j = moments(f)
    fe = feq(rho, j/rho)
    Pi = np.einsum('qk,q...->k...', PAIR_CC, f - fe)
    fneq_reg = (W[:, None, None, None]/(2*CS2**2))*np.einsum('qk,k...->q...', QT*MULT[None, :], Pi)
    return fe + (1 - omega)*fneq_reg


def pressure(rw, ro, G):
    """Shan-Chen 多组分模型的总压力（格子单位）。"""
    return CS2*(rw + ro) + CS2*G*rw*ro


class TwoPhase:
    """油水两相 + 固体格点的小型求解器（与 C++ oilwaterShanChen3D.hh 同一套公式）。

    参数
      shape      : (nx, ny, nz) 网格；一维/二维算例把多余方向设为 1
      tau_w,tau_o: 水、油松弛时间（τ_o = 0.5 + (τ_w−0.5)·M_运动，M_运动 为运动粘度比）
      G          : 跨组分排斥强度
      rho_res    : 另一相的残余密度（C++ 的 RHO_RESIDUAL）
      theta      : 固体接触角（度，从水相一侧量；<90 水湿）；C++ 的 WALL_CONTACT_ANGLE_DEG_TYPEk
      wet_gain   : 接触角外推系数（C++ 的 WETTING_GAIN，标定值 2.4）
      tau_d      : 组分间动量交换松弛时间（C++ 的 DIFFUSION_TAU，默认 1.0；0 = 旧做法）
      wall_model : 'geom'（默认）| 'repel' | 'virtual'（后两种为旧写法，用 s_wall、rho_wall）
      mix_tau    : True=局部混合松弛时间（C++ MIXTURE_TAU=1）；False=旧做法
      body_force : 可选体积力 (gx, gy, gz)，用于分层 Poiseuille 验证
    环境变量（临时切换做对比）：LBM2P_WALL=geom|repel|virtual，LBM2P_TAUD=τ_D，LBM2P_FIXED_TAU=1，
                              LBM2P_TRACE=0（关闭修订八：痕量组分松弛过渡 + 密度地板）
    """

    def __init__(self, shape, tau_w=0.8, tau_o=2.5, G=2.0, rho_res=0.005,
                 theta=90.0, wet_gain=2.4, tau_d=1.0, wall_model='geom',
                 s_wall=0.5, rho_wall=1.0, body_force=(0.0, 0.0, 0.0), u_cap=0.3,
                 mix_tau=True, s_lo=0.05, s_hi=0.95, cap_mode='auto'):
        import os
        self.shape = shape
        self.u_cap = u_cap      # 痕量组分"偏离公共速度"的上限（格子单位）
        self.ow, self.oo = 1.0/tau_w, 1.0/tau_o
        if os.environ.get("LBM2P_FIXED_TAU") == "1":
            mix_tau = False
        self.mix_tau, self.s_lo, self.s_hi = mix_tau, s_lo, s_hi
        self.cap_mode = ('dev' if mix_tau else 'full') if cap_mode == 'auto' else cap_mode
        self.nu_w, self.nu_o = (tau_w - 0.5)/3, (tau_o - 0.5)/3
        self.ncapped = 0
        self.G, self.rho_res = G, rho_res
        self.wall_model = os.environ.get("LBM2P_WALL", wall_model)
        self.tau_d = float(os.environ.get("LBM2P_TAUD", tau_d))
        self.trace = os.environ.get("LBM2P_TRACE", "1") != "0"   # 修订八（C++ 2026-10-27）；LBM2P_TRACE=0 关闭
        self.nrepair = 0
        th = np.radians(theta)
        self.K = wet_gain*np.cos(th)/max(np.sin(th), 1e-6)
        self.s_wall, self.rho_wall = s_wall, rho_wall
        self.g = np.array(body_force, float)
        self.solid = np.zeros(shape, bool)
        self._geom = None
        z = np.zeros((3,) + shape)
        self.fw = feq(np.full(shape, 1.0), z)
        self.fo = feq(np.full(shape, rho_res), z)

    # ----- 初始化工具 -----
    def set_composition(self, water_mask, width=0.0):
        """旧初始化：water_mask 为 True 处放水、其余放油，width>0 时用高斯平滑过渡，总密度不调。"""
        phi = water_mask.astype(float)
        if width > 0:
            from scipy.ndimage import gaussian_filter
            phi = gaussian_filter(phi, sigma=width/2.0, mode='wrap')
        rw = self.rho_res + (1.0 - self.rho_res)*phi
        ro = self.rho_res + (1.0 - self.rho_res)*(1.0 - phi)
        z = np.zeros((3,) + self.shape)
        self.fw, self.fo = feq(rw, z), feq(ro, z)

    def set_composition_pb(self, water_mask, width=2.0):
        """与 C++ 2026-10-26 版一致的初始化：tanh 过渡（半宽 width 格）+ 过渡区总密度按
        "Shan-Chen 压力处处相等"给出（ρ_t + G s(1−s) ρ_t² = P0），避免第 0 步的压力脉冲。"""
        from scipy.ndimage import distance_transform_edt
        m = water_mask.astype(bool)
        d = np.where(m, distance_transform_edt(m), -distance_transform_edt(~m))
        phi = 0.5*(1 + np.tanh((d - 0.5)/max(width, 1e-6))) if width > 0 else m.astype(float)
        r = self.rho_res
        s = r + (1 - 2*r)*phi
        a = self.G*s*(1 - s)
        P0 = 1.0 + self.G*r*(1 - r)
        rt = np.where(a > 1e-12, (np.sqrt(1 + 4*a*P0) - 1)/(2*np.maximum(a, 1e-12)), P0)
        z = np.zeros((3,) + self.shape)
        self.fw, self.fo = feq(s*rt, z), feq((1 - s)*rt, z)

    def set_solid(self, mask):
        self.solid = mask.astype(bool)
        self._geom = None

    # ----- 跨组分力（含壁面润湿）-----
    def _forces(self, rw, ro):
        """返回 Shan-Chen 力 Fw, Fo（不含体积力）。固体格点上的值无意义（会被忽略）。"""
        fl = (~self.solid).astype(float)
        if self.wall_model == 'geom':
            if self._geom is None:
                S = neighbour_sum(self.solid.astype(float))          # Σ t_i c_i（固体邻格）
                Sn = np.sqrt((S**2).sum(0))
                nw = np.where(Sn > 1e-12, -S/np.maximum(Sn, 1e-12), 0.0)   # 法向：固体 → 流体
                nbs = [(i, np.roll(self.solid, shift=tuple(-C[i]), axis=(0, 1, 2))) for i in range(1, Q)]
                self._geom = (S, nw, [(i, m.astype(float)) for i, m in nbs if m.any()])
            S, nw, nbs = self._geom
            rt = rw + ro
            phi = (rw - ro)/rt
            grad = 3.0*(neighbour_sum(phi*fl) + phi[None]*S)        # 固体方向先按镜像 φ(x)
            gt = grad - (grad*nw).sum(0)[None]*nw
            gtm = np.sqrt((gt**2).sum(0))
            Hw = neighbour_sum(rw*fl); Ho = neighbour_sum(ro*fl)
            pm = 1 - 2*self.rho_res
            for i, m in nbs:
                cgt = sum(C[i, d]*gt[d] for d in range(3))
                cn = sum(C[i, d]*nw[d] for d in range(3))
                phis = np.clip(phi + cgt - self.K*gtm*cn, -pm, pm)
                rws, ros = rt*(1 + phis)/2, rt*(1 - phis)/2
                for d in range(3):
                    if C[i, d] != 0:
                        Hw[d] += W[i]*m*rws*C[i, d]
                        Ho[d] += W[i]*m*ros*C[i, d]
            return -self.G*rw*Ho, -self.G*ro*Hw
        # 旧写法：固体格点用固定虚拟密度
        if self.wall_model == 'virtual':
            sw, so = self.s_wall*self.rho_wall, (1.0 - self.s_wall)*self.rho_wall
        else:   # 'repel'
            sw, so = max(0.0, 2*self.s_wall - 1)*self.rho_wall, max(0.0, 1 - 2*self.s_wall)*self.rho_wall
        rw_eff, ro_eff = np.where(self.solid, sw, rw), np.where(self.solid, so, ro)
        return -self.G*rw*neighbour_sum(ro_eff), -self.G*ro*neighbour_sum(rw_eff)

    # ----- 核心：一步时间推进 -----
    def step(self):
        rw, jw = moments(self.fw)
        ro, jo = moments(self.fo)
        if self.trace:     # (0) 修订八(b)：痕量组分密度 < 0.5·ρ_res 时重置为"ρ_res、跟随另一组分速度"的平衡态
            r = self.rho_res
            for (fa, ra), (rb, jb) in (((self.fw, rw), (ro, jo)), ((self.fo, ro), (rw, jw))):
                bad = (ra < 0.5*r) & (rb > r) & ~self.solid
                if bad.any():
                    fa[:] = np.where(bad[None], feq(np.full(self.shape, r), jb/np.maximum(rb, 1e-30)), fa)
                    self.nrepair += int(bad.sum())
            rw, jw = moments(self.fw)
            ro, jo = moments(self.fo)
        rw = np.maximum(rw, self.rho_res)
        ro = np.maximum(ro, self.rho_res)
        # (1)(2) Shan-Chen 跨组分力（含润湿）+ 可选体积力
        Fw, Fo = self._forces(rw, ro)
        Fw = Fw + rw*self.g[:, None, None, None]
        Fo = Fo + ro*self.g[:, None, None, None]
        # (3a) 粘性松弛频率：混合粘度（调和平均）或各自固定
        if self.mix_tau:
            sat = rw/(rw + ro)
            ph = np.clip((sat - self.s_lo)/(self.s_hi - self.s_lo), 0.0, 1.0)
            nu = 1.0/(ph/self.nu_w + (1.0 - ph)/self.nu_o)
            omw = omo = 1.0/(0.5 + 3.0*nu)
        else:
            omw, omo = self.ow, self.oo
        # (3b) 组分间动量交换频率：ω_D=1/τ_D（τ_D=0 时沿用粘性 ω，旧做法）
        exw, exo = (1.0/self.tau_d, 1.0/self.tau_d) if self.tau_d > 0 else (omw, omo)
        # (3c) 修订八(a)：质量分数 < s_lo 的痕量组分，粘性松弛频率逐渐过渡到 ω_D（本体组分不变）
        visc_w, visc_o = omw, omo
        if self.trace and self.mix_tau:
            omt = 1.0/self.tau_d if self.tau_d > 0 else 1.0
            xw = rw/(rw + ro)
            ww = np.minimum(1.0, xw/self.s_lo); wo = np.minimum(1.0, (1.0 - xw)/self.s_lo)
            visc_w, visc_o = ww*omw + (1 - ww)*omt, wo*omo + (1 - wo)*omt
        up = (exw*jw + exo*jo)/(exw*rw + exo*ro)
        # (4) 碰撞后应有的动量 j* = (1−ω_D) j + ω_D ρ u' + F；随后的正则化碰撞保持它
        jw_star = (1 - exw)*jw + exw*rw*up + Fw
        jo_star = (1 - exo)*jo + exo*ro*up + Fo
        # (4') 痕量组分"偏离公共速度"的上限（保护阀）
        for jstar, rho in ((jw_star, rw), (jo_star, ro)):
            if self.cap_mode == 'none':
                continue
            if self.cap_mode == 'full':
                speed = np.sqrt((jstar**2).sum(axis=0))/rho
                scale = np.where(speed > self.u_cap, self.u_cap/np.maximum(speed, 1e-30), 1.0)
                self.ncapped = int(((speed > self.u_cap) & ~self.solid).sum())
                jstar *= scale[None]
            else:
                dev = jstar/rho - up
                speed = np.sqrt((dev**2).sum(axis=0))
                scale = np.where(speed > self.u_cap, self.u_cap/np.maximum(speed, 1e-30), 1.0)
                self.ncapped = int(((speed > self.u_cap) & ~self.solid).sum())
                jstar[:] = rho*(up + dev*scale[None])
        fluid = ~self.solid
        fw = self.fw + (feq(rw, jw_star/rw) - feq(rw, jw/rw))*fluid[None]
        fo = self.fo + (feq(ro, jo_star/ro) - feq(ro, jo/ro))*fluid[None]
        # (5) 碰撞：流体格点正则化 BGK（粘性 ω）；固体格点完全反弹
        fw_c = collide_regularized(fw, visc_w)
        fo_c = collide_regularized(fo, visc_o)
        fw = np.where(self.solid[None], fw[OPP], fw_c)
        fo = np.where(self.solid[None], fo[OPP], fo_c)
        # (6) 迁移
        self.fw, self.fo = stream(fw), stream(fo)

    def run(self, nsteps, check_every=200):
        """推进 nsteps 步；若出现 NaN/Inf 提前返回 False。"""
        for t in range(1, nsteps + 1):
            self.step()
            if t % check_every == 0 and not (np.isfinite(self.fw).all() and np.isfinite(self.fo).all()):
                return False
        return bool(np.isfinite(self.fw).all() and np.isfinite(self.fo).all())

    # ----- 结果提取 -----
    def fields(self):
        """返回 ρ_w, ρ_o, 总速度 u（含半步力修正，固体格点置 0）, 压力 p。"""
        rw, jw = moments(self.fw)
        ro, jo = moments(self.fo)
        Fw, Fo = self._forces(np.maximum(rw, self.rho_res), np.maximum(ro, self.rho_res))
        u = (jw + jo + 0.5*(Fw + Fo))/(rw + ro)
        u = np.where(self.solid[None], 0.0, u)
        return rw, ro, u, pressure(rw, ro, self.G)

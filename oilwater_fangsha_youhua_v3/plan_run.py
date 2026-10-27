#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plan_run.py —— 提交超算前的"可行性计算器"（2026-10-25 新增）
=====================================================================
【它解决什么问题】
  oilwaterDrive 启动后要先建格子、读颗粒，几分钟后才做两项预检：
    ① PRESSURE_PREFLIGHT：岩心两端压差折成格子密度差 Δρ 不能超过 MAX_PREDICTED_LB_DENSITY_DROP；
    ② TIME_BUDGET：MAX_T 必须够注完 TOTAL_PV。
  不通过就退出（你 10 月 25 日的 7 次报错全是这两类）。本脚本用和 C++ 完全相同的公式，
  在登录节点上 1 秒内算完，告诉你：能不能过、V_INF 最大能取多少、要算多少步、大概要多久。

【用法】（在算例目录下）
    python3 plan_run.py                                   # 读 config/defaults.env
    CONFIG_FILE=config/experiments/xxx.env python3 plan_run.py
    python3 plan_run.py --np 64 --mlups-per-core 0.5      # 按 64 核、每核 0.5 MLU/s 估算墙钟时间
    python3 plan_run.py --scan                            # 额外列出 τ_w × 粘度比 的代价表

  --mlups-per-core 怎么定：第一次试算后，run_stdout.log 里每条
      [step=...] t=... PV=... 12.3s 45.6 MLU/s
  的最后一个数是"全部进程合计"的速度，除以 NP 就是每核速度，填进来即可。
  默认 0.5 只是保守猜测（两套 D3Q19 格子 + Shan-Chen + 颗粒耦合）。

【核心公式（与 oilwaterDrive.cpp 一致）】
    Δx = REFERENCE_DIAMETER_UM / N
    ν_lb,w = (TAU_W_TARGET − 0.5)/3 ,   Δt = ν_lb,w · Δx² / ν_w          （ν_w = MU_F/RHO_F × NU_SCALE）
    1 个格子密度单位对应的压力 = ρ_f · c_s² · (Δx/Δt)² ,  c_s² = 1/3
    Δp = V_INF · μ · L / k_min        （达西定律；μ 取油、水中较大者——岩心初始充满油）
    Δρ = Δp / [ρ_f c_s² (Δx/Δt)²]
    注 1 PV 的步数 = 孔隙体积 /(V_INF · 截面积) / Δt
  把 V_INF 取到 Δρ 恰好等于上限，就得到"每 PV 最少步数"：
    步数/PV ≈ 3 · φ · (L²/k) · ν_lb,w · (μ_o/μ_w) / Δρ_max
  —— 注意里面没有 N：加密网格不能减少每 PV 的步数，反而让每步更贵（格点数 ∝ N³）。
"""
from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parameter_config import load_parameters, real, text, integer  # noqa: E402

MD_TO_M2 = 9.869233e-16


def plan(cfg: dict, tau_w: float | None = None, mu_o: float | None = None, v_inf: float | None = None,
         phi: float | None = None) -> dict:
    """按给定参数（可覆盖 τ_w、μ_o、V_INF）算出预检和时间预算的全部数字。"""
    N = integer(cfg, "N")
    ref_d = real(cfg, "REFERENCE_DIAMETER_UM") * 1e-6
    lx, ly, lz = real(cfg, "DOMAIN_LX_M"), real(cfg, "DOMAIN_LY_M"), real(cfg, "DOMAIN_LZ_M")
    core_x0, core_x1 = real(cfg, "CORE_X0_M"), real(cfg, "CORE_X1_M")
    rho_f, mu_f = real(cfg, "RHO_F"), real(cfg, "MU_F")
    rho_o = real(cfg, "RHO_O")
    mu_o = real(cfg, "MU_O") if mu_o is None else mu_o
    v = real(cfg, "V_INF") if v_inf is None else v_inf
    nu_scale = real(cfg, "NU_SCALE")
    tau = real(cfg, "TAU_W_TARGET") if tau_w is None else tau_w
    k_min = real(cfg, "K_MIN_MD") * MD_TO_M2
    phi = 0.5 * (real(cfg, "PHI_MIN") + real(cfg, "PHI_MAX")) if phi is None else phi
    drop_max = real(cfg, "MAX_PREDICTED_LB_DENSITY_DROP")
    total_pv = real(cfg, "TOTAL_PV")
    ramp = real(cfg, "RAMP_TIME")
    equil = real(cfg, "EQUILIBRATION_TIME_S")
    safety = real(cfg, "MIN_TIME_SAFETY_FACTOR")
    max_t = real(cfg, "MAX_T")
    pressure_control = text(cfg, "INLET_CONTROL") == "pressure"
    p_in = real(cfg, "INLET_PRESSURE_PA")
    stop_restart = text(cfg, "DRIVE_PROTOCOL") == "stop_restart"
    hold = real(cfg, "HOLD_TIME_S") + integer(cfg, "RESTART_STAGES") * real(cfg, "RESTART_STAGE_HOLD_S") if stop_restart else 0.0

    nu_w = mu_f / rho_f * nu_scale
    dx = ref_d / N
    if tau > 0.5:
        nu_lb = (tau - 0.5) / 3.0
    else:                                     # TAU_W_TARGET=0：沿用 UMAX（原版做法）
        re = ref_d * v / nu_w
        nu_lb = real(cfg, "UMAX") * N / re
        tau = 0.5 + 3 * nu_lb
    dt = nu_lb * dx * dx / nu_w
    nu_o_lb = nu_lb * (mu_o / rho_o) / (mu_f / rho_f)
    tau_o = 0.5 + 3 * nu_o_lb

    nx, ny, nz = int(round(lx / dx)) + 1, int(round(ly / dx)) + 1, int(round(lz / dx)) + 1
    mx0 = max(int(math.floor(core_x0 / dx + 0.5)), 1)
    mx1 = min(int(math.floor(core_x1 / dx + 0.5)), nx - 2)
    length = max(mx1 - mx0 - 2, 1) * dx                      # 压降测量长度（与 C++ pressureLengthPhys 相同）
    area = (ny - 2) * (nz - 2) * dx * dx
    pore_volume = phi * (mx1 - mx0 + 1) * (ny - 2) * (nz - 2) * dx ** 3

    pa_per_rho = rho_f * (dx / dt) ** 2 / 3.0               # 1 个格子密度单位 = 多少 Pa
    mu_worst = max(mu_f, mu_o)
    factor = max(1.0, real(cfg, "RESTART_MAX_FACTOR")) if stop_restart else 1.0   # 重启最后一级的放大倍数
    if pressure_control:
        dp = factor * p_in
        v_eff = p_in * k_min / (mu_worst * length)          # 恒压：按最粘的相估算流速
    else:
        dp = factor * v * mu_worst * length / k_min
        v_eff = v
    drho = dp / pa_per_rho
    drho_water = v_eff * mu_f * length / k_min / pa_per_rho
    v_max = v_eff * drop_max / drho if drho > 0 else float("inf")

    def drive_time(vel: float) -> float:
        q = vel * area
        target, ramp_vol = total_pv * pore_volume, 0.5 * q * ramp
        t = math.sqrt(2 * target * ramp / q) if target <= ramp_vol else target / q + 0.5 * ramp
        return t + hold

    t_pv = pore_volume / (v_eff * area)
    need_t = safety * (equil + drive_time(v_eff))
    need_t_vmax = safety * (equil + drive_time(v_max))
    return dict(N=N, dx=dx, dt=dt, tau_w=tau, tau_o=tau_o, nu_lb=nu_lb, M=mu_o / mu_f, nx=nx, ny=ny, nz=nz,
                cells=nx * ny * nz, length=length, area=area, pore_volume=pore_volume, phi=phi,
                pa_per_rho=pa_per_rho, dp=dp, drho=drho, drho_water=drho_water, drop_max=drop_max,
                v=v_eff, v_max=v_max, t_pv=t_pv, steps_pv=t_pv / dt, steps_pv_vmax=pore_volume / (v_max * area) / dt,
                total_pv=total_pv, need_t=need_t, need_t_vmax=need_t_vmax, max_t=max_t,
                pressure_control=pressure_control, L2k=length * length / k_min, factor=factor)


def walltime_h(steps: float, cells: int, np_: int, mlups_core: float) -> float:
    """墙钟时间（小时）：总格点更新数 / (核数 × 每核速度)。"""
    return steps * cells / (np_ * mlups_core * 1e6) / 3600.0


def main() -> int:
    ap = argparse.ArgumentParser(description="oilwater_fangsha 提交前可行性估算")
    ap.add_argument("--np", type=int, default=None, help="MPI 进程数（默认读 NP）")
    ap.add_argument("--mlups-per-core", type=float, default=0.5, help="每核速度 MLU/s（试算后用实测值替换）")
    ap.add_argument("--phi", type=float, default=None, help="岩心离散孔隙度（默认 (PHI_MIN+PHI_MAX)/2；N=40 实测 0.3476）")
    ap.add_argument("--scan", action="store_true", help="列出 τ_w × 粘度比 的代价表")
    args = ap.parse_args()

    cfg = load_parameters()
    np_ = args.np or integer(cfg, "NP")
    r = plan(cfg, phi=args.phi)
    ok_drop = r["drho"] <= r["drop_max"]
    ok_time = r["pressure_control"] or r["max_t"] >= r["need_t"]
    mix = integer(cfg, "MIXTURE_TAU")

    print("=" * 78)
    print(f" 配置: {os.environ.get('CONFIG_FILE') or 'config/defaults.env'}")
    print(f" 网格 N={r['N']}  Δx={r['dx']*1e6:.2f} μm  格点 {r['nx']}×{r['ny']}×{r['nz']} = {r['cells']/1e6:.1f} M（×2 套格子）")
    print(f" τ_w={r['tau_w']:.3f}  τ_o={r['tau_o']:.3f}  粘度比 M=μ_o/μ_w={r['M']:.1f}  Δt={r['dt']:.3e} s  MIXTURE_TAU={mix}")
    print(f" 1 个格子密度单位 = {r['pa_per_rho']:.0f} Pa；岩心 L={r['length']*1e3:.3f} mm，L²/k_min={r['L2k']:.3e}")
    print("-" * 78)
    tag = "恒压" if r["pressure_control"] else f"V_INF={r['v']:.3g} m/s"
    print(f" [压差预检] {tag}: 充满水 Δρ≈{r['drho_water']:.3f}；充满{'油' if r['M']>1 else '水'}（最不利）Δρ≈{r['drho']:.3f}"
          f"  上限 {r['drop_max']}  → {'通过' if ok_drop else '不通过'}")
    if not r["pressure_control"]:
        print(f"            满足上限的最大流速 V_INF ≤ {r['v_max']:.3g} m/s")
    else:
        print(f"            满足上限的最大入口压差 INLET_PRESSURE_PA ≤ {r['dp']*r['drop_max']/r['drho']/r['factor']:.4g} Pa")
    print(f" [时间预算] 注 1 PV ≈ {r['t_pv']:.3g} s ≈ {r['steps_pv']:.3g} 步；TOTAL_PV={r['total_pv']} 需 MAX_T ≥ {r['need_t']:.3g} s"
          f"（{r['need_t']/r['dt']:.3g} 步），当前 MAX_T={r['max_t']} → {'通过' if ok_time else '不通过'}")
    if not ok_drop and not r["pressure_control"]:
        print(f"            若把 V_INF 降到 {r['v_max']:.3g}：注 1 PV ≈ {r['steps_pv_vmax']:.3g} 步，"
              f"MAX_T 需 ≥ {r['need_t_vmax']:.3g} s")
    steps_run = (r["need_t"] if ok_drop else r["need_t_vmax"]) / r["dt"]
    h = walltime_h(steps_run, r["cells"], np_, args.mlups_per_core)
    print(f" [墙钟估算] {steps_run:.3g} 步 × {r['cells']/1e6:.1f} M 格点，{np_} 核 × {args.mlups_per_core} MLU/s/核"
          f" ≈ {h:.3g} 小时 ≈ {h/24:.3g} 天")
    if not mix and r["tau_o"] >= 4.5:
        print(" [τ 检查] MIXTURE_TAU=0 且 τ_o ≥ 4.5：C++ 会直接退出（有壁面时发散），请设 MIXTURE_TAU=1。")
    print("=" * 78)

    if args.scan:
        print("\n 代价表（每注 1 PV 的最少步数，V_INF 取到 Δρ 上限；与 N 无关）")
        mus = sorted({1.0e-3 * m for m in (2, 5, 10, 20)} | {real(cfg, "MU_O")})
        head = "  τ_w   " + "".join(f"  M={m/real(cfg,'MU_F'):5.1f}     " for m in mus)
        print(head)
        for tw in (0.7, 0.8):
            row = f"  {tw:.1f}  "
            for m in mus:
                rr = plan(cfg, tau_w=tw, mu_o=m, phi=args.phi)
                row += f"  {rr['steps_pv_vmax']:9.3g}步  "
            print(row)
        print("  （τ_w=0.6 的小液滴、接触角算例在 numpy 验证中发散，不列出；0.6 的数字是 0.7 的一半。\n"
              "   墙钟 ≈ 步数 × 格点数 /(核数 × 每核速度)，格点数 ∝ N³。）")
    return 0 if (ok_drop and ok_time) else 1


if __name__ == "__main__":
    sys.exit(main())

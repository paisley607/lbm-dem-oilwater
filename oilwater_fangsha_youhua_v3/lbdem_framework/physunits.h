/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

/*
 * ============================ 文件说明 ============================
 * physunits.h —— PhysUnits3D：物理单位（SI）与格子单位（lattice units）的换算器。
 *
 * LBM 内部一切都用"格子单位"：格距 Δx = 1、时间步 Δt = 1、参考密度 ρ = 1。
 * LIGGGHTS 用物理单位（m、s、kg）。每次把颗粒位置/速度送进格子、把力送回
 * 颗粒，都要经过本类换算。
 *
 * 换算思路分两层（与 Palabos 的 IncomprFlowParam 一致）：
 *   第 1 层：选一组参考量把物理量无量纲化——
 *            特征长度 l_p（本案例取颗粒直径 2r）、特征速度 u_p、
 *            由此得到特征时间 t_p = l_p/u_p，再加参考密度 ρ_p；
 *   第 2 层：无量纲量再除以无量纲格距 Δx、无量纲时间步 Δt 得到格子单位——
 *            Δx = 1/N（N = 一个 l_p 里有几个格子），Δt = uMax·Δx
 *            （uMax = 特征速度 u_p 对应的格子速度，由用户给定，决定马赫数）。
 *            IncomprFlowParam 还由 Re = l_p·u_p/ν_p 反算出格子粘度与松弛时间 τ。
 *
 * 因此：  物理长度 L → L/l_p/Δx 格子
 *         物理时间 t → t/t_p/Δt 步
 *         物理速度 u → u/u_p·uMax
 * 力的换算因子 forceFactor = ρ_p·(Δx·l_p)⁴/(Δt·t_p)²：
 *   格子单位里"1 个力" = (质量单位)·(长度单位)/(时间单位)²，
 *   质量单位 = ρ_p·(Δx_phys)³，长度单位 = Δx_phys，时间单位 = Δt_phys，
 *   合起来就是 ρ_p·Δx_phys⁴/Δt_phys²，其中 Δx_phys = Δx·l_p，Δt_phys = Δt·t_p。
 * 力矩多乘一个长度 → ρ_p·Δx_phys⁵/Δt_phys²。
 *
 * 另有 0.5 的位置偏移：LIGGGHTS 的域从 0 开始，而 Palabos 第 0 个格点位于
 * 第一个格子的中心 (0.5Δx)，所以位置换算要减/加 0.5（可用 setLbOffset 调整）。
 *
 * 被谁使用：主程序（构造一次）、ibDataExchangeWrappers3D.hh（每步换算）。
 * ================================================================
 */

#ifndef PHYSUNITS_H
#define PHYSUNITS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <cmath>

namespace plb {
  template<typename T>
  class PhysUnits3D {
  public:
    /*
     * 参数（全部物理单位）：
     *   l_p_   特征长度 [m]          u_p_   特征速度 [m/s]
     *   nu_p_  运动粘度 [m²/s]       lx_,ly_,lz_  计算域尺寸 [m]
     *   N      每个特征长度的格子数   uMax   特征速度对应的格子速度（无量纲，通常 ≤ 0.1）
     *   rho_p_ 流体密度 [kg/m³]
     */
    PhysUnits3D(T const l_p_, T const u_p_, T const nu_p_, T const lx_, T const ly_, T const lz_,
                plint const N, T const uMax, T const rho_p_)
      : l_p(l_p_), u_p(u_p_), t_p(l_p_/u_p_), nu_p(nu_p_), lx(lx_), ly(ly_), lz(lz_), rho_p(rho_p_),
        // 雷诺数 Re = l_p·u_p/ν_p；IncomprFlowParam 的后三个参数是用 l_p 无量纲化的域尺寸
        Re(l_p_*u_p_/nu_p_), param(uMax,Re,N,lx_/l_p_,ly_/l_p_,lz/l_p_),
        // forceFactor  = ρ_p · (Δx·l_p)⁴ / (Δt·t_p)²   [kg·m/s² = N]
        forceFactor(pow(param.getDeltaX(),4)/pow(param.getDeltaT(),2) * rho_p * pow(l_p,4)/pow(t_p,2)),
        // torqueFactor = ρ_p · (Δx·l_p)⁵ / (Δt·t_p)²   [N·m]
        torqueFactor(pow(param.getDeltaX(),5)/pow(param.getDeltaT(),2) * rho_p * pow(l_p,5)/pow(t_p,2)),
        lbOffsetX(0.5), lbOffsetY(0.5), lbOffsetZ(0.5)
    {}

    T getRe() const {return Re; };

    // 修改位置偏移（默认各方向 0.5 格），仅在格子与 DEM 域原点不按默认方式对齐时使用
    void setLbOffset(T const lbOffsetX_,T const lbOffsetY_,T const lbOffsetZ_)
    {
      lbOffsetX = lbOffsetX_;
      lbOffsetY = lbOffsetY_;
      lbOffsetZ = lbOffsetZ_;
    }

    // ---------- 物理 → 格子 ----------
    // 长度：L[m] / l_p → 无量纲 → / Δx → 格子数。用于半径等"长度差"，不带偏移
    T getLbLength(T const physL) const { return physL / l_p / param.getDeltaX() ; };
    // 位置：同上再减 0.5，把"域从 0 开始"换成"第 0 格点在 0.5Δx"的 Palabos 约定
    T getLbPosition(T const physL) const { return physL / l_p / param.getDeltaX() - 0.5; };
    T getLbPositionX(T const physL) const { return physL / l_p / param.getDeltaX() - lbOffsetX ; };
    T getLbPositionY(T const physL) const { return physL / l_p / param.getDeltaX() - lbOffsetY ; };
    T getLbPositionZ(T const physL) const { return physL / l_p / param.getDeltaX() - lbOffsetZ ; };
    // 速度：u[m/s] / u_p → 无量纲 → × uMax（= Δx/Δt）→ 格子速度
    T getLbVel(T const physVel) const { return physVel / u_p * param.getLatticeU(); };
    // 时间：t[s] / t_p → 无量纲 → / Δt → 时间步数（实数）
    T getLbTime(T const physTime) const { return physTime / t_p / param.getDeltaT(); };
    // 密度：相对参考密度，ρ_p 本身对应格子密度 1
    T getLbDensity(T const physDensity) const { return physDensity / rho_p; };
    // 加速度 [m/s²]：× t_p²/l_p 无量纲化，再 × Δt²/Δx 进格子（例如重力 g）
    T getLbAccel(T const physAccel) const { return physAccel * (t_p*t_p/l_p)
        * (param.getDeltaT()*param.getDeltaT()/param.getDeltaX()); };
    // 频率/角速度 [1/s 或 rad/s]：× t_p × Δt → 每步的弧度。颗粒角速度 ω 用它换算
    T getLbFreq(T const physFreq) const { return physFreq * t_p * param.getDeltaT(); };
    // 时间 → 整数步数（四舍五入）
    plint getLbSteps(T const physTime) const { return (plint) (physTime/param.getDeltaT()/t_p + 0.5); };
    // 压力 [Pa] → 格子密度偏差 Δρ：
    //   格子里 p_lb = c_s²·Δρ = Δρ/3，压力单位 = 力单位/长度单位²，
    //   故 Δρ = 3·(p − p_offset) / (forceFactor / Δx_phys²)。注意返回的是 ρ−1，不是 ρ
    T getLbRho(T const physPress, T const pOffset = 0) const { return (physPress - pOffset) * 3.
        / getPhysForce(1.) * pow(getPhysLength(1.),2); };

    // ---------- 格子 → 物理（上面各式的逆）----------
    T getPhysLength(T const lbLength) const { return lbLength*param.getDeltaX()*l_p; };
    T getPhysPosition(T const lbLength) const { return (lbLength+0.5)*param.getDeltaX()*l_p; };
    T getPhysPositionX(T const lbLength) const { return (lbLength+lbOffsetX)*param.getDeltaX()*l_p; };
    T getPhysPositionY(T const lbLength) const { return (lbLength+lbOffsetY)*param.getDeltaX()*l_p; };
    T getPhysPositionZ(T const lbLength) const { return (lbLength+lbOffsetZ)*param.getDeltaX()*l_p; };
    T getPhysVel(T const lbVel) const { return lbVel/param.getLatticeU()*u_p; };
    T getPhysTime(T const lbTime) const { return lbTime*param.getDeltaT()*t_p; };
    T getPhysDensity(T const lbDensity) const { return lbDensity * rho_p; };
    T getPhysAccel(T const lbAccel) const { return lbAccel / getLbAccel(1.); };
    // 力/力矩：直接乘换算因子（推导见文件说明）。getForcesFromLattice 用它们把
    // 动量交换得到的格子力送回 LIGGGHTS
    T getPhysForce(T const lbForce) const { return lbForce * forceFactor; };
    T getPhysTorque(T const lbTorque) const { return lbTorque * torqueFactor; };
    // 格子密度 ρ → 物理压力：p = (ρ−1)/3 × 压力单位 + p_offset（与 getLbRho 互逆）
    T getPhysPress(T const lbRho, T pOffset = 0) const { return (lbRho - 1.) / 3. * getPhysForce(1)
        / pow(getPhysLength(1.),2) + pOffset; };


    // 暴露 Palabos 的无量纲参数对象：主程序用它取 omega（松弛频率）、Nx/Ny/Nz 等
    const IncomprFlowParam<T>& getLbParam() const
    {
      return param;
    }

  private:
    T const l_p, u_p, t_p, nu_p;     // 特征长度、速度、时间、运动粘度
    T const lx,ly,lz, rho_p;         // 域尺寸、参考密度
    T const Re;                      // 雷诺数
    IncomprFlowParam<T> param;       // Palabos 无量纲参数（Δx、Δt、uMax、τ…）
    T const forceFactor, torqueFactor;
    T lbOffsetX, lbOffsetY, lbOffsetZ;
  };

};

#endif // PHYSUNITS_H

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
 * ibCompositeDynamics3D.hh —— IBcompositeDynamics 的模板实现，
 * 整个 LBM-DEM 耦合里"流体如何感受到颗粒、颗粒如何受到流体力"的公式都在
 * 本文件的 collide() 里。
 *
 * Noble–Torczynski 部分饱和法（PSM）的碰撞公式（每个离散速度方向 i）：
 *
 *   f_i(x+c_i, t+1) = f_i(x,t) + (1-B)·Ω_i^fluid + B·Ω_i^solid
 *
 *   Ω_i^fluid：普通流体碰撞算子（BGK/Smagorinsky 等，由 baseDynamics 提供）；
 *   Ω_i^solid：固体"非平衡反弹"（non-equilibrium bounce-back）算子：
 *              Ω_i^s = [f_ī − f_ī^eq(ρ,u)] − [f_i − f_i^eq(ρ,u_s)]
 *              其中 ī 为 i 的反方向，u 为格点流体速度，u_s 为颗粒在该点的速度；
 *   B：固体权重函数，B = ε(τ−0.5) / ((1−ε) + (τ−0.5))（LBDEM_USE_WEIGHING 开启时），
 *      ε 为固含率（solid fraction），τ 为松弛时间 = 1/ω。
 *
 * 颗粒受到的水动力（hydrodynamic force）由动量交换得到：
 *   F = − Σ_i c_i · B·Ω_i^s        （格子单位，每格点一份，存入 particleData）
 * 之后由 ibProcessors3D.hh 的 SumForceTorque3D 按颗粒编号汇总并算力矩。
 *
 * 其余函数（序列化、clone、id 注册）是 Palabos 对任何自定义动力学的"例行要求"。
 * ================================================================
 */

#ifndef IB_COMPOSITE_DYNAMICS_HH_LBDEM
#define IB_COMPOSITE_DYNAMICS_HH_LBDEM

#include "ibDef.h"
#include "latticeBoltzmann/dynamicsTemplates.h"

namespace plb {
  // 向 Palabos 的动力学注册表登记本类，名字 "IBcomposite"，换回一个整数 id。
  // 有了 id，Palabos 才能在存盘/并行通信时识别并重建这种动力学；
  // utils.h 也靠这个 id 判断某格点是否是 IB 格点。
  template<typename T, template<typename U> class Descriptor>
  int IBcompositeDynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor, IBcompositeDynamics<T,Descriptor> >("IBcomposite");

  // 三个 static 临时数组的定义（声明见 .h）：全类共享，避免每个格点各存 3×q 个数
  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBcompositeDynamics<T,Descriptor>::fEqSolid =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBcompositeDynamics<T,Descriptor>::fEq =
    Array<T,Descriptor<T>::q>();

  template<typename T, template<typename U> class Descriptor>
  Array<T,Descriptor<T>::q> IBcompositeDynamics<T,Descriptor>::fPre =
    Array<T,Descriptor<T>::q>();

  // 常规构造：把被包装的流体动力学交给父类 CompositeDynamics 保管
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
							 bool automaticPrepareCollision_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_,automaticPrepareCollision_)
  { }

  // 反序列化构造：先建一个空壳（baseDynamics=0），再从字节流里读回全部内容
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(HierarchicUnserializer &unserializer)
    : CompositeDynamics<T,Descriptor>(0,false)
  {
    // pcout << "entering serialize constructor" << std::endl;
    unserialize(unserializer);
  }

  // 拷贝构造：注意 particleData 是普通成员，会被默认拷贝；fEq 等是 static，不需要拷
  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::IBcompositeDynamics(const IBcompositeDynamics &orig)
    : CompositeDynamics<T,Descriptor>(orig)
  { }

  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>::~IBcompositeDynamics() {}

  template<typename T, template<typename U> class Descriptor>
  IBcompositeDynamics<T,Descriptor>* IBcompositeDynamics<T,Descriptor>::clone() const {
    return new IBcompositeDynamics<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  int IBcompositeDynamics<T,Descriptor>::getId() const
  {
    return id;
  }

  /*
   * serialize / unserialize：Palabos 在 MPI 进程间交换"幽灵层"格点、或存盘
   * （checkpoint）时，需要把动力学对象打平成一串数字再恢复。
   * 本类比普通动力学多了 particleData，因此必须把它的 4 个字段也写进去，
   * 顺序为：uPart(d 个) → hydrodynamicForce(d 个) → solidFraction → partId，
   * 最后再交给父类处理 baseDynamics。读回时顺序必须严格一致。
   */
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::serialize(HierarchicSerializer &serializer) const
  {
    for(plint i=0;i<Descriptor<T>::d;i++)
      serializer.addValue(particleData.uPart[i]);

    for(plint i=0;i<Descriptor<T>::d;i++)
      serializer.addValue(particleData.hydrodynamicForce[i]);

    serializer.addValue(particleData.solidFraction);
    serializer.addValue<plint>(particleData.partId);

    CompositeDynamics<T,Descriptor>::serialize(serializer);
  }

  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::unserialize(HierarchicUnserializer &unserializer)
  {
    // 字节流里记录的动力学 id 必须和本类一致，否则说明数据错位
    PLB_PRECONDITION( unserializer.getId() == this->getId() );

    for(plint i=0;i<Descriptor<T>::d;i++)
      unserializer.readValue(particleData.uPart[i]);

    for(plint i=0;i<Descriptor<T>::d;i++)
      unserializer.readValue(particleData.hydrodynamicForce[i]);

    unserializer.readValue(particleData.solidFraction);
    unserializer.readValue<plint>(particleData.partId);

    CompositeDynamics<T,Descriptor>::unserialize(unserializer);
  }

  /*
   * prepareCollision：在真正碰撞之前，把当前格点的全部分布函数 f_i 备份到 fPre。
   * PSM 的固体算子 Ω^s 要用"碰撞前"的 f_i 和 f_ī，而 baseDynamics 的碰撞会
   * 原地改写 cell，所以必须先存一份。只有固含率 > SOLFRAC_MIN（格点被颗粒
   * 触及）时才备份，纯流体格点省掉这次拷贝。
   */
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell)
  {

    if(particleData.solidFraction > SOLFRAC_MIN)
      fPre = cell.getRawPopulations();

  }

  /*
   * defineVelocity：把格点强制设为速度 u、密度 ρ=1 的平衡态。
   * 在 setSpheresOnLattice(initVelFlag=true) 时，被颗粒完全覆盖的格点会用
   * 颗粒速度调用它，避免模拟一开始"静止流体 vs 运动颗粒"产生的受力尖峰。
   * 注意 Palabos 用 rhoBar = ρ-1 表示密度，这里 rhoBar=0 即 ρ=1。
   */
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::defineVelocity(Cell<T,Descriptor>& cell,
                                                         Array<T,Descriptor<T>::d> const& u)
  {
    T const rhoBar = 1.;
    T const uSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(u);
    CompositeDynamics<T,Descriptor>::getBaseDynamics().computeEquilibria(fEq,0.,u,uSqr);
    for(plint i=0;i<Descriptor<T>::q;i++)
      cell[i] = fEq[i];
  }

  /// Implementation of the collision step
  /*
   * collide：PSM 碰撞步，按以下顺序执行
   *   (1) 备份碰撞前分布 fPre；
   *   (2) 清零本格点的水动力累加器；
   *   (3) 若格点不是"完全固体"，先做一次普通流体碰撞（得到 Ω^fluid）；
   *   (4) 若格点是"纯流体"（ε < SOLFRAC_MIN），到此结束；
   *   (5) 否则计算两套平衡分布：fEq（流体速度）和 fEqSolid（颗粒速度）；
   *   (6) 按 ε 落在哪个区间，分别用"完全反弹"或"B 加权混合"公式更新分布，
   *       并同时累加动量交换力。
   */
  template<typename T, template<typename U> class Descriptor>
  void IBcompositeDynamics<T,Descriptor>::collide(Cell<T,Descriptor>& cell,
                                                  BlockStatistics& statistics)
  {

    // (1) 备份 f_pre（仅当格点被颗粒触及时才真正拷贝）
    prepareCollision(cell);

    // (2) 本格点对颗粒的力每步重新累加，先清零
    particleData.hydrodynamicForce.resetToZero();

    // (3) ε <= SOLFRAC_MAX：格点内还有流体，执行被包装的普通碰撞。
    //     碰撞后 cell[i] = fPre[i] + Ω_i^fluid，后面用 (cell - fPre) 反推出 Ω^fluid。
    //     ε > SOLFRAC_MAX（完全在颗粒内部）则跳过，流体碰撞对它没有意义。
    if(particleData.solidFraction <= SOLFRAC_MAX)
      CompositeDynamics<T,Descriptor>::collide(cell,statistics);
      // this->getBaseDynamics().collide(cell,statistics);

    // (4) 纯流体格点：普通碰撞已经完成，无需任何 PSM 修正
    if(particleData.solidFraction < SOLFRAC_MIN)
      return;

    // (5a) 由碰撞后的分布求宏观量：rhoBar = ρ-1，j = ρu（动量密度）
    T rhoBar(0.);
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell,rhoBar,j);

    // (5b) fEq = f^eq(ρ, u)：以流体速度计算的二阶（Ma^2）平衡分布
    T const jSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(j);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                         j, jSqr, fEq );

    // (5c) fEqSolid = f^eq(ρ, u_s)：以颗粒速度计算的平衡分布。
    //      bgk_ma2_equilibria 接收的是动量 j=ρu 而不是 u，所以要把 uPart 乘以 ρ=(1+rhoBar)。
    Array<T,Descriptor<T>::d> uPart = particleData.uPart*(1.+rhoBar);
    T const uPartSqr = VectorTemplateImpl<T,Descriptor<T>::d>::normSqr(uPart);
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar),
                                                         uPart, uPartSqr, fEqSolid );


    if(particleData.solidFraction > SOLFRAC_MAX){
      /*
       * (6a) 完全固体格点：相当于 B=1，f_i^new = f_i^pre + Ω_i^s，
       *      即只做"非平衡反弹"，不含任何流体碰撞成分。
       *
       * 静止方向 i=0（c_0=0，反方向是它自己）单独处理：Ω_0^s = f_0^eq(u_s) − f_0^eq(u)。
       * ★★ 2026-10-26 修正：上游 LBDEMcoupling 原写法是 fEq[0] − fEqSolid[0]（符号反了），
       *    每步凭空增减质量 Δρ = ρ(u_s² − u²)（把 Σ_i Ω_i^s 加起来即可验证：正确写法恰为 0）。
       *    单相流 u 很小，这个误差看不出来；两相流里"痕量组分"（油区里只有 0.005 的水）的
       *    速度 j/ρ 可以到 0.1–1，颗粒内部格点每步就丢掉与自身密度同量级的质量，
       *    几步之内密度变成 −1（云端复现：第 8 步颗粒内部 ρ_w≈−1.0，随后整场 NaN，
       *    骨架颗粒受力 1.8e9 N）。改成正确符号后两相质量严格守恒。
       */
      T const bb0 = fEqSolid[0] - fEq[0];
      cell[0] = fPre[0] + bb0;

      // Palabos 的 D3Q19 中，方向 i (1..q/2) 的反方向是 i+q/2，因此一次循环
      // 同时处理一对相反方向 (iPop, iOpp)。
      for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
        plint const iOpp = iPop + Descriptor<T>::q/2;
        // bb    = Ω_iPop^s = [f_ī − f_ī^eq(u)] − [f_i − f_i^eq(u_s)]
        // bbOpp = Ω_iOpp^s，把 i 和 ī 互换即可
        T const bb = ((fPre[iOpp]-fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));
        T const bbOpp = ((fPre[iPop]-fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));

        cell[iPop] = fPre[iPop] + bb;
        cell[iOpp] = fPre[iOpp] + bbOpp;

        // 动量交换：流体在这对方向上获得的动量为 c_i·bb + c_ī·bbOpp = c_i·(bb−bbOpp)
        // （因 c_ī = −c_i）。颗粒受力 = 流体动量变化的负值，故用 "-="。
        for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
          particleData.hydrodynamicForce[iDim]
            -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
      }
    } else { /* particleData.solidFraction < SOLFRAC_MAX */
      /*
       * (6b) 部分覆盖格点（SOLFRAC_MIN < ε <= SOLFRAC_MAX）：流体碰撞与固体反弹
       *      按权重 B 混合。
       */
#ifdef LBDEM_USE_WEIGHING
      // B = ε(τ−0.5) / ((1−ε) + (τ−0.5))，τ = 1/ω 为松弛时间。
      // ooo 即 (τ−0.5)，它正比于运动粘度 ν = (τ−0.5)/3。
      // 这种加权使 ε→1 时 B→1、ε→0 时 B→0，中间段随粘度调整，减弱结果对 τ 的敏感性。
      T const ooo = 1./CompositeDynamics<T,Descriptor>::getBaseDynamics().getOmega() - 0.5;
      T const B = particleData.solidFraction*ooo/((1.- particleData.solidFraction) + ooo);
#else
      // 不加权：直接用固含率作为固体权重
      T const B = particleData.solidFraction;
#endif
      T const oneMinB = 1. - B;

      // 静止方向：f_0^new = f_0^pre + (1−B)·Ω_0^fluid + B·Ω_0^s，Ω_0^s = f_0^eq(u_s) − f_0^eq(u)，
      // 其中 Ω_0^fluid = cell[0] − fPre[0]（普通碰撞已在步骤 (3) 把结果写进 cell）。
      // ★★ 2026-10-26 修正符号（同上，原写法每步造成质量误差 B·ρ(u_s² − u²)）。
      T const bb0 = fEqSolid[0] - fEq[0];
      cell[0] = fPre[0] + oneMinB*(cell[0] - fPre[0]) + B*bb0;

      for(plint iPop=1;iPop<=Descriptor<T>::q/2;iPop++){
        plint const iOpp = iPop + Descriptor<T>::q/2;
        // 这里的 bb 已经乘上了 B，即 bb = B·Ω_i^s
        T const bb = B*((fPre[iOpp] - fEq[iOpp]) - (fPre[iPop]-fEqSolid[iPop]));
        T const bbOpp = B*((fPre[iPop] - fEq[iPop]) - (fPre[iOpp]-fEqSolid[iOpp]));


        // 完整的 PSM 更新式：f_i^new = f_i^pre + (1−B)·Ω_i^fluid + B·Ω_i^s
        cell[iPop] = fPre[iPop] + oneMinB*(cell[iPop] - fPre[iPop]) + bb;
        cell[iOpp] = fPre[iOpp] + oneMinB*(cell[iOpp] - fPre[iOpp]) + bbOpp;

        // 动量交换力 F = −Σ_i c_i·B·Ω_i^s，与 (6a) 同理
        for(plint iDim=0;iDim<Descriptor<T>::d;iDim++)
          particleData.hydrodynamicForce[iDim]
            -= Descriptor<T>::c[iPop][iDim]*(bb-bbOpp);
      }
    } /* if solidFraction > SOLFRAC_MAX */

  } /* IBcompositeDynamics::collide() */

}; /* namespace plb */

#endif /* IB_COMPOSITE_DYNAMICS_HH_LBDEM */

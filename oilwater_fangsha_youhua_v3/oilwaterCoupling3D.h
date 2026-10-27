/*
 * oilwaterCoupling3D.h
 * ====================================================================================
 * 两相 LBM → DEM（LIGGGHTS）的力交换
 *
 * 【给初学者的说明】
 *   颗粒受到的流体力由三部分组成，本文件把它们加在一起写入 LIGGGHTS 的耦合 fix
 *   （couple/lb/onetoone），DEM 再用这些力积分颗粒运动：
 *     (1) 水格子的 PSM 动量交换力（压力 + 粘性，由 IBcompositeDynamics 逐格累计）；
 *     (2) 油格子的 PSM 动量交换力；
 *     (3) 毛细反作用力：Shan-Chen 排斥在"流体—颗粒"之间产生的力对颗粒的反作用，
 *         由 OilWaterShanChenFunctional3D 累加在 TwoPhaseForceBuffers 里。
 *   原版只有 (1)+(2)，缺 (3)，也就是审计报告维度 2 指出的"毛细力没有传到颗粒上"。
 *
 * 【为什么不能直接调用框架自带的 getForcesFromLattice 两次】
 *   框架的 getForcesFromLattice 每次都会先把耦合 fix 里的力清零再写入，
 *   对水、油各调一次的话后一次会把前一次覆盖掉。所以这里统一：清零一次 →
 *   三部分累加 → 转物理单位 → 写 fix → 一次 MPI 归约。
 * ====================================================================================
 */

#ifndef OILWATER_COUPLING_3D_H
#define OILWATER_COUPLING_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "plb_ib.h"
#include "oilwaterShanChen3D.h"   // TwoPhaseForceBuffers

namespace plb {

/// 把水、油两套格子的 PSM 力 + 毛细反作用力相加，写入 DEM 耦合 fix。
/// 返回本进程写入的颗粒数（含幽灵颗粒，调试用）。
template <typename T, template <typename U> class Descriptor>
plint getTwoPhaseForcesFromLattice(MultiBlockLattice3D<T, Descriptor>& latticeW,
                                   MultiBlockLattice3D<T, Descriptor>& latticeO,
                                   LiggghtsCouplingWrapper& wrapper,
                                   PhysUnits3D<T> const& units,
                                   TwoPhaseForceBuffers<T> const* capillary,
                                   T& maxCapillaryForcePhys);

}  // namespace plb

#include "oilwaterCoupling3D.hh"

#endif  // OILWATER_COUPLING_3D_H

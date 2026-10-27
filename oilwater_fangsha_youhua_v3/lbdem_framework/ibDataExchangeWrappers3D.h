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
 * some wrapper functions for data exchange between LIGGGHTS and Palabos
 */

/*
 * ============================ 文件说明 ============================
 * ibDataExchangeWrappers3D.h —— 主程序直接调用的两个"高层耦合接口"的声明。
 *
 * 一个耦合步的典型流程（见 oilwaterDrive.cpp 主循环）：
 *   1. setSpheresOnLattice(lattice, wrapper, units, ...)
 *        把 LIGGGHTS 当前所有颗粒的位置/速度映射到格子（写 particleData）；
 *   2. lattice.collideAndStream()  若干次
 *        流体演化，碰撞步同时算出每个格点对颗粒的力；
 *   3. getForcesFromLattice(lattice, wrapper, units)
 *        把格点力汇总成每个颗粒的合力/合力矩，换成物理单位，交给 LIGGGHTS；
 *   4. wrapper.run(n)
 *        DEM 推进若干步，颗粒在流体力+接触力作用下运动；回到 1。
 *
 * 本文件只做声明；实现在 ibDataExchangeWrappers3D.hh。
 * 底层工作由 ibProcessors3D 的 SetSingleSphere3D / SumForceTorque3D 完成，
 * 单位换算由 physunits.h 的 PhysUnits3D 完成。
 *
 * setSpheresOnLattice 有三个重载：
 *   (a) 带 excludeType：跳过某些 LIGGGHTS 颗粒类型（例如壁面颗粒）；
 *   (b) 不带 excludeType：最常用的简化版；
 *   (c) 本案例新增，带 couplingDomain：只在给定的全局盒子内做固相映射。
 * ================================================================
 */

#ifndef IBDATAEXCHANGEWRAPPERS_H_LBDEM
#define IBDATAEXCHANGEWRAPPERS_H_LBDEM

#include "ibProcessors3D.h"

namespace plb {
  /*
   * if initVelFlag is set true, the populations at the grid nodes covered with
   * particles are initialized with the particle velocity. This is useful to avoid
   * force spikes when a particle already has a certain velocity upon beginning of
   * the simulation. During normal simulation, it should be set to false
   */
  /*
   * initVelFlag：为 true 时，被颗粒完全覆盖的格点的分布函数会被初始化为颗粒
   * 速度对应的平衡态。只在模拟开始、颗粒已有初速度时用一次；正常耦合循环中
   * 必须传 false，否则每步都会强行重置颗粒内部的流体状态。
   *
   * excludeType：LIGGGHTS 颗粒类型列表，属于这些类型的颗粒不映射到格子。
   */
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           std::vector<plint> &excludeType,
                           bool initVelFlag);

  // 简化版：不排除任何颗粒类型，映射范围不限制
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag);

  /*
   * 带耦合域（couplingDomain）的重载：固含率只映射到耦合域内的格子。
   * 用于把出口压力边界（及其相邻层）排除在 IBM 固相映射之外，避免球体
   * 跨骑 Zou-He 等压边界时破坏边界动量平衡。couplingDomain 为全局格子
   * 坐标；传 0 表示不限制（等价于 4 参版本）。
   */
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag,
                           Box3D const* couplingDomain);

  /*
   * getForcesFromLattice：汇总格点上的水动力到颗粒，换算成物理单位后写入
   * LIGGGHTS 的 fix couple/lb/onetoone，并做 MPI 归约使每个颗粒的力完整。
   * 调用前必须已经执行过至少一次 collideAndStream（否则 hydrodynamicForce 全为 0）。
   */
  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units);

}; /* namespace plb */

#include "ibDataExchangeWrappers3D.hh"

#endif /* IBDATAEXCHANGEWRAPPERS_H_LBDEM */

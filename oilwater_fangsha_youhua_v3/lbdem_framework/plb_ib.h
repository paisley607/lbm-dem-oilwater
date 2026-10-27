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
 * aggregate header for ib in plb
 */

/*
 * ============================ 文件说明 ============================
 * plb_ib.h —— 耦合框架的"总入口"头文件（aggregate header）。
 *
 * 主程序（oilwaterDrive.cpp）只需 #include "plb_ib.h" 就能拿到整个
 * LBM-DEM 耦合框架。它本身不含代码，只按依赖顺序汇总以下模块：
 *
 *   ibCompositeDynamics3D.h     —— PSM/IBM 碰撞算子（流体侧核心）
 *   ibProcessors3D.h            —— 固含率映射 SetSingleSphere3D、受力汇总 SumForceTorque3D
 *   ibDataWritingFunctionals3D.h—— 把格点上的颗粒信息导出成场，供 VTK 后处理
 *   ibDataExchangeWrappers3D.h  —— 高层封装：setSpheresOnLattice / getForcesFromLattice
 *   latticeDecomposition.h      —— 让 Palabos 的并行分块与 LIGGGHTS 的域分解一致
 *   physunits.h                 —— 物理单位 <-> 格子单位换算
 *
 * 其中 ibProcessors3D.h 又会包含 liggghtsCouplingWrapper.h（LIGGGHTS 封装）、
 * ibDef.h（全局宏）和 utils.h（从格点取颗粒数据的工具函数）。
 * periodicPressureFunctionals3D.h 不在此列，需要时由主程序单独包含（本案例未用）。
 * ================================================================
 */

#ifndef PLB_IB_H
#define PLB_IB_H

#include "ibCompositeDynamics3D.h"
#include "ibProcessors3D.h"
#include "ibDataWritingFunctionals3D.h"
#include "ibDataExchangeWrappers3D.h"
#include "latticeDecomposition.h"
#include "physunits.h"

#endif // PLB_IB_H

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
 * helper class to create lattice decomposition from LIGGGHTS parallelization
 */

/*
 * ============================ 文件说明 ============================
 * latticeDecomposition.h —— LatticeDecomposition 类的声明：
 * 让 Palabos 的格子并行分块（block decomposition）与 LIGGGHTS 的空间域分解
 * （domain decomposition）完全重合。
 *
 * 为什么必须一致？
 *   耦合时，每个 MPI 进程要把"自己看得见的颗粒"（本地+幽灵）映射到
 *   "自己拥有的格子块"上，再把格子块上的力收回给颗粒。只有当某进程的
 *   LIGGGHTS 子域与它的 Palabos 格子块在空间上是同一块区域时，这种
 *   "本进程内部直接读写、不额外通信"的做法才成立。否则颗粒在 A 进程、
 *   它覆盖的格子却在 B 进程，就得做复杂的跨进程数据交换。
 *
 * 做法：LIGGGHTS 启动后已经决定了 procgrid（各方向进程数）和 xsplit/ysplit/zsplit
 * （各方向的切分比例）；本类读出这些信息，按同样比例把 nx×ny×nz 的格子切块，
 * 并用 ExplicitThreadAttribution 把每个块显式指派给对应的 MPI rank。
 *
 * 主程序用法：
 *   LatticeDecomposition lDec(nx,ny,nz,wrapper.lmp);
 *   MultiBlockLattice3D<T,DESCRIPTOR> lattice(
 *       MultiBlockManagement3D(lDec.getBlockDistribution(), lDec.getThreadAttribution(), envelopeWidth),
 *       defaultMultiBlockPolicy3D().getBlockCommunicator(),
 *       defaultMultiBlockPolicy3D().getCombinedStatistics(),
 *       defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(), dynamics);
 *
 * 实现见 latticeDecomposition.cpp（非模板，需单独编译链接）。
 * ================================================================
 */

#ifndef LATTICE_DECOMPOSITION_H
#define LATTICE_DECOMPOSITION_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "lammps.h"

namespace plb{

class LatticeDecomposition {
public:
  // nx_,ny_,nz_：全局格子尺寸；lmp_：已初始化（已读入 LIGGGHTS 脚本、域已建立）的 LAMMPS 对象
  LatticeDecomposition(plb::plint nx_, plb::plint ny_, plb::plint nz_, LAMMPS_NS::LAMMPS *lmp_);

  ~LatticeDecomposition();

  // 返回分块结构的副本（按值）——描述每个块覆盖哪个 Box3D
  plb::SparseBlockStructure3D getBlockDistribution();
  // 返回新分配的线程指派对象——描述每个块归哪个 MPI rank；调用者（MultiBlockManagement3D）接管所有权
  plb::ExplicitThreadAttribution* getThreadAttribution();
private:
  plb::plint nx,ny,nz;
  LAMMPS_NS::LAMMPS &lmp;
  plb::plint npx,npy,npz;                    // LIGGGHTS 各方向的进程数
  std::vector<plb::plint> xVal, yVal, zVal;  // 各方向的切分点（格子坐标），长度 np+1
  plb::SparseBlockStructure3D *blockStructure;
  plb::ExplicitThreadAttribution *threadAttribution;
};

}

#endif /* LATTICE_DECOMPOSITION_H */

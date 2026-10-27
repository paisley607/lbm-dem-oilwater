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
 * latticeDecomposition.cpp —— LatticeDecomposition 的实现：
 * 把 LIGGGHTS 的域分解"翻译"成 Palabos 的格子分块。
 *
 * LIGGGHTS（LAMMPS）侧的并行信息都在 lmp->comm 里：
 *   procgrid[3]      —— x/y/z 方向各有几个进程（npx×npy×npz = 总进程数）；
 *   xsplit/ysplit/zsplit —— 长度 np+1 的数组，给出各方向切分面的相对位置
 *                       （0 = 域左端，1 = 域右端），例如 2 进程均分时 xsplit = {0, 0.5, 1}；
 *   grid2proc[i][j][k] —— 三维进程网格坐标 → MPI rank 的映射。
 *
 * 翻译规则：
 *   切分面相对位置 × 格子尺寸，四舍五入 → 格子坐标的切分点 xVal/yVal/zVal；
 *   第 (i,j,k) 个块 = [xVal[i], xVal[i+1]-1] × [yVal[j], yVal[j+1]-1] × [zVal[k], zVal[k+1]-1]，
 *   指派给 rank grid2proc[i][j][k]。
 *
 * 这样得到的块恰好与各进程的 LIGGGHTS 子域重合，且每个进程恰好一个块——
 * ibDataExchangeWrappers3D.hh 里 getBlocks()[0] 的用法正依赖这一点。
 *
 * 前提：LIGGGHTS 的模拟盒（region/create_box）必须与 Palabos 的格子覆盖同一
 * 物理区域，且 PhysUnits3D 的换算让两者坐标原点一致。
 * ================================================================
 */

#include "latticeDecomposition.h"

// LAMMPS / LIGGGHTS includes
#include "comm.h"   // 提供 lmp.comm->procgrid / xsplit / grid2proc

namespace plb {

LatticeDecomposition::LatticeDecomposition(plb::plint nx_, plb::plint ny_, plb::plint nz_, LAMMPS_NS::LAMMPS *lmp_)
  : nx(nx_),ny(ny_),nz(nz_),lmp(*lmp_),
    blockStructure(0),threadAttribution(0)
{
  // 第 1 步：读出 LIGGGHTS 各方向的进程数
  npx = lmp.comm->procgrid[0];
  npy = lmp.comm->procgrid[1];
  npz = lmp.comm->procgrid[2];

  // 第 2 步：把相对切分位置（0~1）换成格子坐标的整数切分点。
  // 共 np+1 个点：xVal[0]=0，xVal[npx]=nx。round 避免浮点截断导致漏掉或重复一层格子。
  for(plint i=0;i<=npx;i++)
    xVal.push_back(round( lmp.comm->xsplit[i]*(double)nx ));
  for(plint i=0;i<=npy;i++)
    yVal.push_back(round( lmp.comm->ysplit[i]*(double)ny ));
  for(plint i=0;i<=npz;i++)
    zVal.push_back(round( lmp.comm->zsplit[i]*(double)nz ));

  // 第 3 步：建立覆盖整个格子 [0,nx-1]×[0,ny-1]×[0,nz-1] 的稀疏块结构
  //（Palabos 的 Box3D 是闭区间，所以上界要减 1）
  blockStructure = new SparseBlockStructure3D(Box3D(xVal[0], xVal.back()-1,
                                                    yVal[0], yVal.back()-1,
                                                    zVal[0], zVal.back()-1) );

  threadAttribution = new ExplicitThreadAttribution();

  // 第 4 步：逐个生成块并指派进程。块的遍历顺序 (iX,iY,iZ) 与 LIGGGHTS 的
  // 进程网格坐标一一对应，grid2proc 给出该网格位置上的 MPI rank。
  for (plint iX=0; iX<xVal.size()-1; ++iX) {
    for (plint iY=0; iY<yVal.size()-1; ++iY) {
      for (plint iZ=0; iZ<zVal.size()-1; ++iZ) {
        plint id = blockStructure->nextIncrementalId();
        blockStructure->addBlock (Box3D( xVal[iX], xVal[iX+1]-1, yVal[iY],
                                         yVal[iY+1]-1, zVal[iZ], zVal[iZ+1]-1 ),
                                  id);
        threadAttribution->addBlock(id,(plint)lmp.comm->grid2proc[iX][iY][iZ]);
      }
    }
  }

}

LatticeDecomposition::~LatticeDecomposition()
{
  if(blockStructure != 0) delete blockStructure;
  if(threadAttribution != 0) delete threadAttribution;
}


// 按值返回副本：MultiBlockManagement3D 会自己保存一份
SparseBlockStructure3D LatticeDecomposition::getBlockDistribution()
{
  return SparseBlockStructure3D(*blockStructure);
}
// 返回堆上的新副本：MultiBlockManagement3D 接管并负责 delete，
// 本对象析构时只释放自己持有的那份
ExplicitThreadAttribution* LatticeDecomposition::getThreadAttribution()
{
  return new ExplicitThreadAttribution(*threadAttribution);
}

}

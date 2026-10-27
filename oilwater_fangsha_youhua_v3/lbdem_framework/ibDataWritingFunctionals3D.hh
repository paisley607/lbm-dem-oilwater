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
 * ibDataWritingFunctionals3D.hh —— 两个后处理 functional 的实现。
 *
 * 两者结构完全相同：三重循环遍历 domain 内每个格点 →
 *   用 getParticleDataFromCell 找到该格点的 particleData →
 *   按 which 选出要的字段 → 写进输出场的同一位置。
 * 找不到 particleData（非 IB 格点）时写 0。
 *
 * getTypeOfModification 的约定：modified[0] 对应 lattice（未修改），
 * modified[1] 对应输出场（修改了静态变量，即场的数值），Palabos 据此决定
 * 是否同步输出场的幽灵层。
 *
 * 本文件只在写输出时被调用，对耦合计算结果没有任何影响。
 * ================================================================
 */


#ifndef IBDATAWRITINGFUNCTIONALS3D_HH_LBDEM
#define IBDATAWRITINGFUNCTIONALS3D_HH_LBDEM

#include "utils.h"

namespace plb {
  /*
     implementation of getScalarQuantityFromDynamicsFunctional
  */

  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>::
  process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
          ScalarField3D<T2>& data)
  {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          T1 val = 0; // 默认 0：非 IB 格点或无颗粒覆盖

          Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);

          // 剥开动力学找到颗粒信息（见 utils.h）
          IBdynamicsParticleData<T1,Descriptor>* particleData =
            getParticleDataFromCell<T1,Descriptor>(cell);

          if(particleData){
            switch(which){
            case SolidFraction:   // 固含率 ε，0~1
              val = (T2) particleData->solidFraction;
              break;
            case ParticleId:      // 颗粒编号，0 表示无颗粒；可视化时可按 id 上色
              val = (T2) particleData->partId;
              break;

            }
          }
          // 输出场与格子使用相同的局部坐标，直接写同一位置
          data.get(iX,iY,iZ) = val;
        }
      }
    }
  }
  template<typename T1, template<typename U> class Descriptor, typename T2>
  GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>*
  GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>::clone() const
  {
    return new GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>(*this);
  }

  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>::
  getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;          // 格子只读
    modified[1] = modif::staticVariables;  // 输出场被写入
  }

  /*
     implementation of getVectorQuantityFromDynamicsFunctional
  */

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  void GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>::
  process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
          TensorField3D<T2,nDim>& data)
  {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

          Array<T2,nDim> val;
          val.resetToZero(); // 默认零向量

          Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);

          IBdynamicsParticleData<T1,Descriptor>* particleData =
            getParticleDataFromCell<T1,Descriptor>(cell);

          if(particleData){
            switch(which){
            case ParticleVelocity:   // 颗粒在该格点的速度（含转动贡献），格子单位
              val = particleData->uPart;
              break;
            case HydrodynamicForce:  // 该格点最近一次碰撞算出的动量交换力，格子单位
              val = particleData->hydrodynamicForce;
              break;

            }
          }
          data.get(iX,iY,iZ) = val;
        }
      }
    }
  }

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>*
  GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>::clone() const
  {
    return new GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>(*this);
  }

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  void GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>::
  getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
  }


}; /* namespace plb */

#endif /* IBDATAWRITINGFUNCTIONALS3D_HH_LBDEM */

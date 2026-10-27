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
 * functionals used to obtain particle data that can be used for post processing
 */

/*
 * ============================ 文件说明 ============================
 * ibDataWritingFunctionals3D.h —— 后处理用：把藏在每个格点动力学对象里的
 * 颗粒信息（particleData）"抄"到普通的标量场/张量场里，以便写 VTK 可视化。
 *
 * 背景：Palabos 的 VtkImageOutput3D 只会输出 ScalarField3D / TensorField3D，
 * 而固含率、颗粒 id、颗粒速度、水动力都存在 IBcompositeDynamics::particleData
 * 里，VTK 写不出来。本文件提供两个 data processor 做这个转换：
 *
 *   GetScalarQuantityFromDynamicsFunctional —— 取标量：SolidFraction 或 ParticleId
 *   GetVectorQuantityFromDynamicsFunctional —— 取向量：ParticleVelocity 或 HydrodynamicForce
 *
 * 用法示例（主程序）：
 *   MultiScalarField3D<T> sf(lattice);
 *   applyProcessingFunctional(
 *       new GetScalarQuantityFromDynamicsFunctional<T,DESCRIPTOR,T>(SolidFraction),
 *       lattice.getBoundingBox(), lattice, sf);
 *   vtkOut.writeData<float>(sf, "solidFraction", 1.);
 *
 * 不参与耦合计算，只影响输出。实现见 .hh，依赖 utils.h 的 getParticleDataFromCell。
 * ================================================================
 */

#ifndef IBDATAWRITINGFUNCTIONALS3D_H_LBDEM
#define IBDATAWRITINGFUNCTIONALS3D_H_LBDEM

namespace plb {

  // 可导出的量的枚举：传给构造函数选择要抄哪一个字段
  enum IBscalarQuantity { SolidFraction, ParticleId };
  enum IBvectorQuantity { ParticleVelocity, HydrodynamicForce };

  // _LS = Lattice + Scalar field：process 同时拿到格子和一个标量场
  // T1：格子的浮点类型；T2：输出标量场的类型（可以不同，例如写 float 省空间）
  template<typename T1, template<typename U> class Descriptor, typename T2>
  struct GetScalarQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2> {
  public:

    IBscalarQuantity which;
    GetScalarQuantityFromDynamicsFunctional(IBscalarQuantity const which_)
      : which(which_) {}

    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         ScalarField3D<T2>& data);

    GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

  };

  // _LT = Lattice + Tensor field：nDim 为向量维数（3）
  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  struct GetVectorQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LT<T1,Descriptor,T2,nDim> {
  public:

    IBvectorQuantity which;
    GetVectorQuantityFromDynamicsFunctional(IBvectorQuantity const which_)
      : which(which_) {}

    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         TensorField3D<T2,nDim>& data);

    GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  };

}; /* namespace plb */

#include "ibDataWritingFunctionals3D.hh"

#endif /* IBDATAWRITINGFUNCTIONALS3D_H_LBDEM */

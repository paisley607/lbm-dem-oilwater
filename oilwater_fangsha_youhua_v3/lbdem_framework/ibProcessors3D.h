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
 * ibProcessors3D.h —— 三个 Palabos "data processor"（数据处理器）的声明，
 * 它们是 DEM 颗粒与 LBM 格子之间双向数据搬运的具体执行者。
 *
 * Palabos 里对格子做逐点操作的标准方式是写一个 BoxProcessingFunctional3D_L
 * 子类，实现 process(domain, lattice)，再用 applyProcessingFunctional 施加到
 * 某个 Box3D 区域上；Palabos 负责把区域切到各 MPI 进程的本地块上执行。
 *
 *   SetSingleSphere3D   —— "颗粒 → 格子"：给定一个球（位置 x、速度 v、角速度 ω、
 *                          半径 r、编号 id），对包围盒内每个格点算固含率
 *                          （solid fraction）并写入 particleData。
 *   AttributeFunctional —— 空操作，只用来触发一次 Palabos 的幽灵层同步，
 *                          让各进程边界上的 particleData 保持一致。
 *   SumForceTorque3D    —— "格子 → 颗粒"：遍历本地格点，把碰撞步累加的
 *                          hydrodynamicForce 按颗粒编号求和成合力，并用
 *                          (格点−质心)×力 累加力矩。
 *   ParticleData        —— 只是两个 typedef，给"颗粒位置数组"起个短名字。
 *
 * 调用者：ibDataExchangeWrappers3D.hh 的 setSpheresOnLattice / getForcesFromLattice。
 * 实现：ibProcessors3D.hh。
 * ================================================================
 */
#ifndef IB_PROCESSORS_3D_H
#define IB_PROCESSORS_3D_H

#include "liggghtsCouplingWrapper.h"
#include "physunits.h"
#include "ibDynamicsParticleData.h"

namespace plb{
  // 便于书写的类型别名：ParticleDataArrayVector = 每个颗粒一个 3 维向量的数组
  template<typename T>
  struct ParticleData {
    typedef typename std::vector<Array<T,3> > ParticleDataArrayVector;
    typedef typename std::vector<T> ParticleDataScalarVector;
  };

  /*
   * SetSingleSphere3D：把一个球形颗粒"画"到格子上。
   * 所有几何量（x、v、omega、com、r）都已换算为格子单位（见 physunits.h）。
   * com 是质心坐标——对单球就是球心 x；对多球刚体（multisphere）是整体质心，
   * 用来计算 ω×r 时的 r。
   */
  template<typename T, template<typename U> class Descriptor>
  struct SetSingleSphere3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SetSingleSphere3D(T *x_, T *v_, T *omega_, T *com_, T r_, int id_, bool initVelFlag_ = false)
      : x(x_),v(v_),omega(omega_),com(com_),r(r_),id(id_), initVelFlag(initVelFlag_) {}
    SetSingleSphere3D(const SetSingleSphere3D &orig)
      : x(orig.x),v(orig.v),omega(orig.omega),com(orig.com),
        r(orig.r),id(orig.id),initVelFlag(orig.initVelFlag) {}

    // 颗粒的包围盒（格子坐标）：球心 ± (r+2)。多出的 2 个格子是安全余量，
    // 保证球面附近做 5x5x5 子采样时相关格点都被包进来。
    // 调用者只在这个盒子上施加本 functional，而不是整个格子，以省计算量。
    Box3D getBoundingBox() const { return Box3D(x[0]-r-2,x[0]+r+2,
                                                x[1]-r-2,x[1]+r+2,
                                                x[2]-r-2,x[2]+r+2);}

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    SetSingleSphere3D<T,Descriptor>* clone() const;
  private:
    T *x,*v,*omega,*com,r; // com: center of mass
    int id;
    bool initVelFlag;      // 为真时，被完全覆盖的格点用颗粒速度初始化分布函数
    // 计算以 (dx,dy,dz) 为"格点相对球心偏移"的单位立方体被半径 r 的球占据的体积比
    T calcSolidFraction(T const dx_, T const dy_, T const dz_, T const r_);
    // 把固含率、颗粒速度（含转动）、颗粒 id 写入格点的 particleData
    void setValues(IBdynamicsParticleData<T,Descriptor> &p,
                   T const sf, T const dx, T const dy, T const dz);
    // 清空格点的 particleData（表示无颗粒覆盖）
    void setToZero(IBdynamicsParticleData<T,Descriptor> &p);
  };

  /*
   * this functional does nothing, it is just used to get
   * all the static variables right after looping over
   * the spheres
   */
  /*
   * AttributeFunctional：process() 什么都不做；关键在 getTypeOfModification 返回
   * modif::dynamicVariables——这会让 Palabos 在施加它之后，对各进程块之间的
   * 幽灵层（envelope/ghost layer）做一次同步，把动力学对象里的 particleData
   * 复制到邻居进程。SetSingleSphere3D 自己声明"不修改"（避免每放一个球就同步一次），
   * 所以所有球放完后统一用它同步一次。
   */
  template<typename T, template<typename U> class Descriptor>
  struct AttributeFunctional : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice) {}
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {modified[0] = modif::dynamicVariables;}

    AttributeFunctional<T,Descriptor>* clone() const
    { return new AttributeFunctional<T,Descriptor>(*this);}
  };

  /*
   * SumForceTorque3D：把格点上的水动力汇总到颗粒。
   *   x       —— 所有本地+幽灵颗粒的位置（格子单位），按 LIGGGHTS 本地索引排列；
   *   force/torque —— 输出数组，长度 3*nPart，按本地索引排列，本 functional 向其累加；
   *   wrapper —— 用来通过 atom->map(tag) 把颗粒全局编号换成本地索引。
   */
  template<typename T, template<typename U> class Descriptor>
  struct SumForceTorque3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SumForceTorque3D(typename ParticleData<T>::ParticleDataArrayVector &x,
                     T *force_, T *torque_, LiggghtsCouplingWrapper &wrapper_);
    SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig);

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);

    SumForceTorque3D<T,Descriptor>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

  private:
    typename ParticleData<T>::ParticleDataArrayVector &x;
    T *force,*torque;
    LiggghtsCouplingWrapper &wrapper;
    void addForce(plint const partId, plint const coord, T const value);
    void addTorque(plint const partId, plint const coord, T const value);
  };




}; // plb

#include "ibProcessors3D.hh"

#endif // IB_PROCESSORS_3D_H

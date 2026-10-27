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
 * periodicPressureFunctionals3D.h —— "周期边界 + 压差驱动"的工具类声明。
 *
 * 【本案例（oilwaterDrive.cpp）未使用本文件】：本案例的出入口采用 Zou-He
 * 速度/压力边界而非周期边界，所以这里的类都不会被实例化。保留它是为了
 * 与上游 LBDEMcoupling 完整一致，下面只做简要说明。
 *
 * 问题背景：想模拟"沿 x 方向无限长的管道/多孔介质中的压差驱动流"时，常用
 * 周期边界让流体从出口回到入口；但普通周期边界两端压力相同，没有驱动力。
 * 本文件提供两种在周期边界两侧人为制造压力（密度）跳变的方法：
 *   ZhangPeriodicPressureFunctional3D —— Zhang & Kwok (2006)：按 ρ_target/ρ_avg
 *                                        直接缩放穿越边界的分布函数；
 *   KimPeriodicPressureFunctional3D   —— Kim & Pitsch (2007)：用"减去旧平衡态、
 *                                        加上目标密度的新平衡态"修正分布函数；
 *   PeriodicPressureManager           —— 封装 Kim 方法的完整用法（preColl/postColl）。
 * 另有两个初始化辅助类：
 *   PressureGradient                      —— 用线性密度梯度初始化整个域；
 *   PoiseuilleProfileAndPressureGradient  —— 线性密度梯度 + 抛物线速度剖面初始化。
 * BoxProcessingFunctional3D_LST 是 Kim 方法需要的"格子+标量场+张量场"三元 functional 基类。
 * ================================================================
 */

#ifndef PERIODIC_PRESSURE_FUNCTIONALS_3D
#define PERIODIC_PRESSURE_FUNCTIONALS_3D

namespace plb {
  /*
   * implements a periodic pressure gradient across a periodic boundary
   * following
   * Zhang and Kwok, Phys.Rev.E 73, 047702 (2006)
   */

  // Zhang 方法：构造时筛出"沿 dimension 方向、分量等于 direction（+1 或 -1）"的
  // 那些离散速度（即穿越边界进入域内的方向），process 里只缩放这些分布函数
  template<typename T, template<typename U> class Descriptor>
  class ZhangPeriodicPressureFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    ZhangPeriodicPressureFunctional3D(T rhoTarget_, T rhoAvg_, plint dimension_, plint direction_);
    void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual ZhangPeriodicPressureFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const;

  private:
    typedef std::vector<plint> IndexVec;
    IndexVec rescalePop;        // 需要缩放的离散速度方向的编号列表
    T rhoTarget,rhoAvg;         // 目标密度、边界层当前平均密度
    plint dimension, direction; // 边界法向（0=x,1=y,2=z）与进入方向（±1）

  };

  /**
   * this derives directly from BoxProcessingFunctional
   * it creates a coupling between three objects: a lattice,
   * a scalar field and a tensor field.
   */

  // Palabos 自带 _L、_LS、_LT 等 functional 基类，但没有"格子+标量+张量"三合一的，
  // 这里自己派生一个：processGenericBlocks 把三个 AtomicBlock 转型后转给 process
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  class BoxProcessingFunctional3D_LST : public BoxProcessingFunctional3D {
  public:
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor> &lattice,
                         ScalarField3D<T2> &scalar, TensorField3D<T3,nDim> &tensor) =0;
    void processGenericBlocks (Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);

  };

  /*
   * implements a periodic pressure gradient across a periodic boundary
   * following
   * Kim and Pitsch, Phys. Fluids 19, 108101(2007)
   */

  // Kim 方法：需要对侧边界的密度场 rho 和速度场 u（由 PeriodicPressureManager 事先
  // 交换好），对每个进入方向 i：f_i += f_i^eq(ρ_target_cell, u) − f_i^eq(ρ_per, u)
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  class KimPeriodicPressureFunctional3D : public BoxProcessingFunctional3D_LST<T1,Descriptor,T2,T3,nDim> {
  public:
    KimPeriodicPressureFunctional3D(T1 rhoTarget_, T1 rhoAvg_, plint dimension_, plint direction_);
    void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                 ScalarField3D<T2> &rho, TensorField3D<T3,nDim> &u);
    virtual KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const;

  private:
    typedef std::vector<plint> IndexVec;
    IndexVec rescalePop;
    T1 rhoTarget,rhoAvg;
    plint dimension, direction;
    // 用给定 (ρ, u) 计算第 iEq 个方向的平衡分布（内部把 u 换成动量 j=ρu）
    T1 computeEquilibrium(Cell<T1,Descriptor>& cell,
                          plint iEq,
                          T2 rho,
                          Array<T3,Descriptor<T1>::d> u);

  };

  /*
   * convenience wrapper class that contains instantiation and application
   * of KimPeriodicPressureFunctional3D
   * usage:
   * ppm.preColl(lattice);
   * lattice.collideAndStream();
   * ppm.postColl(lattice);
   * dimension: 0 = x, 1 = y, 2 = z
   * in/out direction: +1 if boundary
   */

  // 管理器：preColl 在碰撞前读取入口/出口层的 ρ、u 并互换（入口用出口的、出口用入口的），
  // postColl 在碰撞-迁移后对两层各施加一次 Kim functional，使入口密度趋向 rhoIn、出口趋向 rhoOut
  template<typename T, template<typename U> class Descriptor>
  class PeriodicPressureManager {
  public:
    PeriodicPressureManager(MultiBlockLattice3D<T,Descriptor> &lattice,
                            T rhoIn_, T rhoOut_, Box3D inlet_, Box3D outlet_,
                            plint dimension_, plint inDirection_, plint outDirection_);

    void preColl(MultiBlockLattice3D<T,Descriptor> &lattice);
    void postColl(MultiBlockLattice3D<T,Descriptor> &lattice);
  private:
    T rhoIn, rhoOut, rhoAvgIn, rhoAvgOut;   // 目标密度与当前层平均密度
    Box3D inlet,outlet,tmp;                 // 入口层、出口层、交换用的临时层
    plint dimension, inDirection, outDirection;
    MultiScalarField3D<T> rho;              // 与格子同尺寸的密度场缓存
    MultiTensorField3D<T,Descriptor<T>::d> u; // 与格子同尺寸的速度场缓存
  };

  /*
   * class PressureGradient
   * helper class to initialize a domain with a constant pressure gradient
   * usage:
   * initializeAtEquilibrium( lattice, lattice.getBoundingBox(),
   *                          PressureGradient<T>(pHi,pLo, nGrid, dimension) );
   * nGrid: grid points in direction of gradient
   * dimension: 0 = x, 1 = y, 2 = z
   */

  // 初始化函子：密度从 pHi 线性降到 pLo（这里 pHi/pLo 实际是格子密度），速度为 0
  template<typename T>
  class PressureGradient {
  public:
    PressureGradient(T pHi_, T pLo_, plint n_, plint dimension_);
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const;
  private:
    T pHi, pLo;
    plint n,dimension;
  };

  // 初始化函子：线性密度梯度 + 沿 dimension 方向的泊肃叶（抛物线）速度剖面
  template<typename T>
  class PoiseuilleProfileAndPressureGradient {
  public:
    PoiseuilleProfileAndPressureGradient(T pHi_, T pLo_,  T uMax,
                                         plint nx_, plint ny_, plint nz_, plint dimension_);
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const;
  private:
    T pHi,pLo,uMax;
    plint nx,ny,nz,dimension;
    T poiseuilleVelocity(T i, T n) const; // 归一化抛物线 1 − ((i − n/2)/(n/2))²
  };

}; /* namespace plb */

#include "periodicPressureFunctionals3D.hh"

#endif /* PERIODIC_PRESSURE_FUNCTIONALS_3D */

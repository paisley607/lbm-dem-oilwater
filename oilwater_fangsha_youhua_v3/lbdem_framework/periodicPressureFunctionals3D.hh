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
 * periodicPressureFunctionals3D.hh —— 周期压差边界各类的模板实现。
 *
 * 【本案例未使用】，详见 .h 的说明；此处只对每段实现给出简要注释，
 * 便于以后需要"周期边界 + 压差驱动"时参考。
 *
 * 两种方法的共同点：只修改"从边界进入计算域"的那几个离散速度方向
 * （rescalePop），使边界层的密度被拉向目标值 rhoTarget，从而在周期边界
 * 两侧形成压力跳变 Δp = (ρ_in − ρ_out)/3。
 *   Zhang：f_i ← (f_i + t_i)·(ρ_target/ρ_avg) − t_i   （t_i 为格子权重，
 *          因 Palabos 存的是 f_i − t_i，先加回去再缩放）
 *   Kim  ：f_i ← f_i − f_i^eq(ρ_per,u_per) + f_i^eq(ρ_per + ρ_target − ρ_avg, u_per)
 *          其中 ρ_per、u_per 取自周期对侧边界层。
 *
 * PeriodicPressureManager 负责 Kim 方法所需的"对侧数据交换"：
 *   preColl：算出入口/出口层的 ρ、u，借助 tmp 层把两者对调，并记录各层平均密度；
 *   postColl：对入口层施加 (rhoIn, rhoAvgOut)，对出口层施加 (rhoOut, rhoAvgIn)。
 * ================================================================
 */

namespace plb {
  /*
   * implementation of BoxProcessingFunctional3D_LST
   */
  // 把 Palabos 传来的三个通用块按约定顺序（格子、标量场、张量场）转型后调用 process
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  void BoxProcessingFunctional3D_LST<T1,Descriptor,T2,T3,nDim>::processGenericBlocks
  (Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks)
  {
    process( domain,
             dynamic_cast<BlockLattice3D<T1,Descriptor>&>(*atomicBlocks[0]),
             dynamic_cast<ScalarField3D<T2>&>(*atomicBlocks[1]),
             dynamic_cast<TensorField3D<T3,nDim>&>(*atomicBlocks[2])
             );
  }

  /*
   * implementation of KimPeriodicPressureFunctional3D
   */

  // 构造：挑出 c_i 在 dimension 方向分量等于 direction 的所有离散速度
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>::KimPeriodicPressureFunctional3D
  (T1 rhoTarget_, T1 rhoAvg_, plint dimension_, plint direction_)
    : rhoTarget(rhoTarget_), rhoAvg(rhoAvg_),
      dimension(dimension_), direction(direction_)
  {
    for(plint iPop=0;iPop<Descriptor<T1>::q;iPop++){
      if(Descriptor<T1>::c[iPop][dimension] == direction){
        rescalePop.push_back(iPop);
      }
    }
  }
  // process：对边界层每个格点，用对侧的 (ρ_per, u_per) 做"减旧平衡、加新平衡"修正
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  void KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>::process
  (Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
   ScalarField3D<T2> &rho, TensorField3D<T3,nDim> &u)
  {
    // std::cout << global::mpi().getRank() << " "
    //           << domain.x0 << " " << domain.x1 << " | "
    //           << domain.y0 << " " << domain.y1 << " | "
    //           << domain.z0 << " " << domain.z1 << std::endl;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);
          // rho/u 场里存的已经是"对侧边界层"的值（由 PeriodicPressureManager::preColl 交换）
          T2 rhoPer = rho.get(iX,iY,iZ);
          Array<T3,3> uPer = u.get(iX,iY,iZ);

          // 目标密度 = 对侧局部密度 + (本侧目标平均 − 对侧实际平均)，保留横向密度分布形状
          T2 rhoTargetCell = rhoPer + (rhoTarget - rhoAvg);


          for(IndexVec::iterator it = rescalePop.begin();
              it != rescalePop.end(); it++){
            plint ii = *it;
            cell[ii] -= computeEquilibrium(cell,ii,rhoPer,uPer);
            cell[ii] += computeEquilibrium(cell,ii,rhoTargetCell,uPer);
          }

        }
      }
    }
  }

  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>*
  KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>::clone() const
  {
    return new KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>(*this);
  }

  // 只修改格子（modified[0]），rho/u 两个场只读
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  void KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>::getTypeOfModification
  (std::vector<plb::modif::ModifT>& modified) const
  {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;//modif::staticVariables;
    modified[2] = modif::nothing;//modif::staticVariables;
  }

  // Palabos 的 Cell::computeEquilibrium 要的是动量 j=ρu 和 |j|²，这里原地把 u 换成 j
  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  T1 KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>::computeEquilibrium(Cell<T1,Descriptor>& cell,
                                                                                   plint iEq,
                                                                                   T2 rho,
                                                                                   Array<T3,Descriptor<T1>::d> u)
  {
    T1 jSqr = 0;
    for(plint i=0;i<Descriptor<T1>::d;i++){
      u[i] *= rho;
      jSqr += u[i]*u[i]; // u now holds j (momentum) but no additional variable needs to be declared
    }
    return cell.computeEquilibrium(iEq,rho,u,jSqr);
  }


  /*
   * implementation of ZhangPeriodicPressureFunctional
   */


  // 构造：同 Kim，挑出进入方向
  template<typename T, template<typename U> class Descriptor>
  ZhangPeriodicPressureFunctional3D<T,Descriptor>::ZhangPeriodicPressureFunctional3D
  (T rhoTarget_, T rhoAvg_, plint dimension_, plint direction_)
    : rhoTarget(rhoTarget_), rhoAvg(rhoAvg_),
      dimension(dimension_), direction(direction_)
  {
    for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
      if(Descriptor<T>::c[iPop][dimension] == direction){
        rescalePop.push_back(iPop);
      }
    }
    // for(plint i=0;i<rescalePop.size();i++)
    //   pcout << rescalePop[i] << " ";
    // pcout << std::endl;
  }

  // process：按 ρ_target/ρ_avg 的比例缩放进入方向的分布函数
  template<typename T, template<typename U> class Descriptor>
  void ZhangPeriodicPressureFunctional3D<T,Descriptor>::process
  (Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    T rescaleFactor = rhoTarget/rhoAvg;
    // pcout << dimension << " " << direction << " | "
    // 	  << rhoTarget << " " << rhoAvg << " " << rescaleFactor << std::endl;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

	  Cell<T,Descriptor> &cell = lattice.get(iX,iY,iZ);

          for(plint i=0;i<rescalePop.size();i++){
	    plint ii=rescalePop[i];
            // Palabos 内部存的是 f_i − t_i（t_i 为格子权重），先加回 t_i 得真实 f_i，
            // 缩放后再减掉 t_i 存回去
	    T fTmp = cell[ii] + Descriptor<T>::t[ii];
	    fTmp *= rescaleFactor;
	    cell[ii] = fTmp - Descriptor<T>::t[ii];
          }

        }
      }
    }
  }

  template<typename T, template<typename U> class Descriptor>
  ZhangPeriodicPressureFunctional3D<T,Descriptor>*
  ZhangPeriodicPressureFunctional3D<T,Descriptor>::clone() const
  {
    return new ZhangPeriodicPressureFunctional3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void ZhangPeriodicPressureFunctional3D<T,Descriptor>::getTypeOfModification
  (std::vector<plb::modif::ModifT>& modified) const
  {
    modified[0] = modif::staticVariables;
  }

  /*
   * implementation of PeriodicPressureManager
   */


  // 构造：tmp 是一个与入口层同形状、但位于坐标 1 处的"中转层"，preColl 用它来交换两层数据
  template<typename T, template<typename U> class Descriptor>
  PeriodicPressureManager<T,Descriptor>::PeriodicPressureManager(MultiBlockLattice3D<T,Descriptor> &lattice,
                                                                 T rhoIn_, T rhoOut_, Box3D inlet_, Box3D outlet_,
                                                                 plint dimension_, plint inDirection_,
                                                                 plint outDirection_)
    : rhoIn(rhoIn_), rhoOut(rhoOut_), rhoAvgIn(0.), rhoAvgOut(0.),
      inlet(inlet_), outlet(outlet_), tmp(inlet_),
      dimension(dimension_), inDirection(inDirection_), outDirection(outDirection_),
    rho(lattice),u(lattice)
  {
    switch(dimension){
    case 0:
      tmp = Box3D(1,1,inlet.y0,inlet.y1,inlet.z0,inlet.z1);
      break;
    case 1:
      tmp = Box3D(inlet.x0,inlet.x1,1,1,inlet.z0,inlet.z1);
      break;
    case 2:
      tmp = Box3D(inlet.x0,inlet.x1,inlet.y0,inlet.y1,1,1);
      break;
    }
  }

  // preColl：计算入口/出口层的 ρ、u，然后三步交换（inlet→tmp, outlet→inlet, tmp→outlet），
  // 使 rho/u 场中入口位置存的是出口数据、出口位置存的是入口数据；再记录两层平均密度
  template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureManager<T,Descriptor>::preColl(MultiBlockLattice3D<T,Descriptor> &lattice)
  {
    computeDensity(lattice,rho,inlet); computeDensity(lattice,rho,outlet);
    computeVelocity(lattice,u,inlet); computeVelocity(lattice,u,outlet);

    copy(rho,inlet,rho,tmp); copy(rho,outlet,rho,inlet); copy(rho,tmp,rho,outlet);
    copy(u,inlet,u,tmp); copy(u,outlet,u,inlet); copy(u,tmp,u,outlet);

    rhoAvgIn = computeAverageDensity(lattice,inlet);
    rhoAvgOut = computeAverageDensity(lattice,outlet);
  }

  // postColl：对入口层施加 Kim 修正（目标 rhoIn，参考对侧平均 rhoAvgOut），出口层同理
  template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureManager<T,Descriptor>::postColl(MultiBlockLattice3D<T,Descriptor> &lattice)
  {
    std::vector<MultiBlock3D*> vecIn, vecOut;
    vecIn.push_back(&lattice); vecIn.push_back(&rho); vecIn.push_back(&u);
    vecOut.push_back(&lattice); vecOut.push_back(&rho); vecOut.push_back(&u);

    applyProcessingFunctional
      (new KimPeriodicPressureFunctional3D<T,Descriptor,T,T,Descriptor<T>::d>(rhoIn,
                                                                              rhoAvgOut,dimension,inDirection),
       inlet,vecIn);
    applyProcessingFunctional
      (new KimPeriodicPressureFunctional3D<T,Descriptor,T,T,Descriptor<T>::d>(rhoOut,
                                                                              rhoAvgIn,dimension,outDirection),
       outlet,vecOut);

  }

  /*
   * implementation of PressureGradient
   */
  template<typename T>
  PressureGradient<T>::PressureGradient(T pHi_, T pLo_, plint n_, plint dimension_)
    : pHi(pHi_), pLo(pLo_), n(n_), dimension(dimension_)
  { }
  // 沿 dimension 方向从 pHi 线性过渡到 pLo，速度置零；
  // 配合 initializeAtEquilibrium 使用，每个格点调用一次
  template<typename T>
  void PressureGradient<T>::operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
  {
    velocity.resetToZero();
    switch(dimension){
    case 0:
      density = pHi - (pHi-pLo)*(T)iX/(T)(n-1);
      break;
    case 1:
      density = pHi - (pHi-pLo)*(T)iY/(T)(n-1);
      break;
    case 2:
      density = pHi - (pHi-pLo)*(T)iZ/(T)(n-1);
      break;
    }
  }


  template<typename T>
  PoiseuilleProfileAndPressureGradient<T>::PoiseuilleProfileAndPressureGradient(T pHi_, T pLo_,  T uMax_,
                                                                                plint nx_, plint ny_, plint nz_,
                                                                                plint dimension_)
    : pHi(pHi_), pLo(pLo_), uMax(uMax_),
      nx(nx_), ny(ny_), nz(nz_), dimension(dimension_)
  { }
  // 线性密度梯度 + 流向速度 = uMax × 两个横向的抛物线剖面之积（矩形管道近似）
  template<typename T>
  void PoiseuilleProfileAndPressureGradient<T>::operator()
    (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
  {
    velocity.resetToZero();
    switch(dimension){
    case 0:
      density = pHi - (pHi-pLo)*(T)iX/(T)(nx-1);
      velocity[0] = uMax*poiseuilleVelocity(iY,ny)*poiseuilleVelocity(iZ,nz);
      velocity[1] = 0.;
      velocity[2] = 0.;
      break;
    case 1:
      density = pHi - (pHi-pLo)*(T)iY/(T)(ny-1);
      velocity[0] = 0.;
      velocity[1] = uMax*poiseuilleVelocity(iX,nx)*poiseuilleVelocity(iZ,nz);
      velocity[2] = 0.;
      break;
    case 2:
      density = pHi - (pHi-pLo)*(T)iZ/(T)(nz-1);
      velocity[0] = 0.;
      velocity[1] = 0.;
      velocity[2] = uMax*poiseuilleVelocity(iX,nx)*poiseuilleVelocity(iY,ny);
      break;
    }
  }
  // 归一化抛物线：中心为 1，两壁（i=0 与 i=n）为 0
  template<typename T>
  T PoiseuilleProfileAndPressureGradient<T>::poiseuilleVelocity(T i, T n) const
  {
    T nHalf = n/2.;
    T x = (i-nHalf)/(nHalf);
    return (1.-x*x);
  }
};

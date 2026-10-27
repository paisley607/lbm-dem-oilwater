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
 * ibProcessors3D.hh —— SetSingleSphere3D 与 SumForceTorque3D 的模板实现。
 *
 * 这是"颗粒 ↔ 格子"双向映射的算法核心：
 *
 * 颗粒 → 格子（SetSingleSphere3D）
 *   对包围盒内每个格点：
 *     1. 用 calcSolidFraction 计算固含率 ε：先用"内切/外接球"快速判断
 *        完全在内（ε=1）或完全在外（ε=0），只有球面附近的格点才做
 *        5×5×5=125 个子采样点的计数；
 *     2. 根据 (新 ε 是否>0, 旧 ε 是否>0) 组成的 decFlag 决定是写入、清零还是保留，
 *        以正确处理"颗粒离开格点"和"两个颗粒共享一个格点"的情况；
 *     3. 写入的颗粒速度 = 平动速度 v + ω × (格点 − 质心)。
 *
 * 格子 → 颗粒（SumForceTorque3D）
 *   对本地块内每个格点：若 partId>0，把 hydrodynamicForce 累加到该颗粒的合力，
 *   并用 r × F（r = 格点 − 颗粒中心）累加力矩；r 要做最小镜像处理以兼容周期边界。
 *
 * 所有量均为格子单位；换算回物理单位在 ibDataExchangeWrappers3D.hh 里完成。
 * ================================================================
 */

#include "ibDef.h"
#include "ibCompositeDynamics3D.h"
#include "ibDynamicsParticleData.h"
#include "utils.h"

// 以下是 LIGGGHTS 的头文件：需要用 atom->map() 做"全局编号 → 本地索引"换算
#include "lammps.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

namespace plb{

  /*
   * implementation of SetSingleSphere3D
   */

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    // modifies stuff, but we don't want updates running
    // after each particle set... should be modif::dataStructure
    // modified[0] = modif::dataStructure;
    //
    // 这里故意"撒谎"说没修改任何东西：否则每放一个颗粒 Palabos 都会做一次
    // 进程间幽灵层通信，颗粒多时代价巨大。所有颗粒放完后，由 AttributeFunctional
    // 统一触发一次同步（见 setSpheresOnLattice）。
    modified[0] = modif::nothing;
  }

  template<typename T, template<typename U> class Descriptor>
  SetSingleSphere3D<T,Descriptor>* SetSingleSphere3D<T,Descriptor>::clone() const
  {
    return new SetSingleSphere3D<T,Descriptor>(*this);
  }

  /*
   * process：对 domain（已经是包围盒 ∩ 本地块）内每个格点更新 particleData。
   * 注意 Palabos 传进来的 iX,iY,iZ 是"块内局部坐标"，要加上 lattice.getLocation()
   * 才是全局格子坐标，而颗粒位置 x[] 是全局格子坐标。
   */
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);

          // 找到该格点动力学里的颗粒信息；非 IB 格点（如边界节点）直接跳过
          IBdynamicsParticleData<T,Descriptor>* particleData =
            getParticleDataFromCell<T,Descriptor>(cell);

          if(!particleData) continue;

          // this one actually helps, believe it or not
          // 提示 CPU 预取这块内存（第二个参数 1 = 准备写入），纯性能优化
          __builtin_prefetch(particleData,1);

          // 局部坐标 → 全局格子坐标
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);


          // 格点相对球心的偏移（用于算固含率）
          T const dx = xGlobal - x[0];
          T const dy = yGlobal - x[1];
          T const dz = zGlobal - x[2];

          // 格点相对质心的偏移（用于算 ω×r；单球时 com == x）
          T const dx_com = xGlobal - com[0];
          T const dy_com = yGlobal - com[1];
          T const dz_com = zGlobal - com[2];

          // 本颗粒在该格点的固含率
          T const sf = calcSolidFraction(dx,dy,dz,r);

          // 该格点上一次（可能是别的颗粒）留下的固含率和颗粒 id
          T const sf_old = particleData->solidFraction;
          int const id_old = (int) particleData->partId;

          /*
           * decFlag 的二进制编码：bit0 = 新固含率>0，bit1 = 旧固含率>0。
           *   0：新旧都为 0 —— 该格点本来就没颗粒，本颗粒也没覆盖它，清零保险；
           *   1：新>0，旧=0 —— 本颗粒新覆盖了一个空格点，直接写入；
           *   2：新=0，旧>0 —— 本颗粒不覆盖它，但格点里有记录：
           *        若记录的是本颗粒（id_old==id），说明颗粒已经移走，清零；
           *        若记录的是别的颗粒，不能动，否则会抹掉别人的信息；
           *   3：新旧都>0 —— 两种情况：
           *        记录就是本颗粒（id_old==id）：更新为最新值；
           *        记录是别的颗粒：两颗粒共享此格点，谁的固含率大听谁的
           *        （sf > sf_old 才覆盖），避免接触点处来回抖动。
           */
          plint const decFlag = (sf > SOLFRAC_MIN) + 2*(sf_old > SOLFRAC_MIN);

          switch(decFlag){
          case 0: // sf == 0 && sf_old == 0
            setToZero(*particleData);
            break; // do nothing
          case 1: // sf > 0 && sf_old == 0
            setValues(*particleData,sf,dx_com,dy_com,dz_com);
            break;
          case 2: // sf == 0 && sf_old > 0
            if( id_old == id ) // then particle has left this cell
              setToZero(*particleData);
            break; // else do nothing
          case 3: // sf > 0 && sf_old > 0
            if( sf > sf_old || id_old == id )
              setValues(*particleData,sf,dx_com,dy_com,dz_com);
            break; // else do nothing
          }
          // if desired, initialize interior of sphere with sphere velocity
          // 仅在初始化阶段使用：把颗粒内部格点的分布函数设成颗粒速度的平衡态，
          // 避免第一步出现巨大的受力尖峰（见 ibDataExchangeWrappers3D.h 的说明）
          if(initVelFlag && sf > SOLFRAC_MAX)
            cell.defineVelocity(particleData->uPart);

        }
      }
    }

  }

  /*
   * calcSolidFraction：计算以格点为中心、边长为 1 的立方体被半径 r 的球占据的体积比。
   * (dx_,dy_,dz_) 为格点相对球心的偏移。
   * 精确的球-立方体相交体积没有简单解析式，这里用"子采样计数"近似：
   * 把单位立方体均分为 5×5×5 = 125 个小格，数一数有多少个小格中心落在球内。
   */
  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFraction(T const dx_, T const dy_, T const dz_, T const r_)
  {
    // return calcSolidFractionRec(dx_,dy_,dz_,r_);

    static plint const slicesPerDim = 5;                       // 每个方向切 5 片
    static T const sliceWidth = 1./((T)slicesPerDim);          // 每片宽 0.2
    static T const fraction = 1./((T)(slicesPerDim*slicesPerDim*slicesPerDim)); // 每个子格占 1/125

    // should be sqrt(3.)/2.
    // add a little to avoid roundoff errors
    // 单位立方体中心到顶点的距离是 sqrt(3)/2 ≈ 0.866；用 sqrt(3.1)/2 略放大一点作保险
    static const T sqrt3half = (T) sqrt(3.1)/2.;

    // 格点到球心距离的平方（避免开方）
    T const dist = dx_*dx_ + dy_*dy_ + dz_*dz_;

    // 快速判定 1：格点中心离球心比 r + 0.866 还远 → 整个立方体都在球外，ε=0
    T const r_p = r_ + sqrt3half;
    if (dist > r_p*r_p) return 0;

    // 快速判定 2：格点中心离球心比 r − 0.866 还近 → 整个立方体都在球内，ε=1
    T const r_m = r_ - sqrt3half;
    if (dist < r_m*r_m) return 1;

    // 只有球面附近的"过渡层"格点才走到这里做 125 点采样
    T const r_sq = r_*r_;
    T dx_sq[slicesPerDim],dy_sq[slicesPerDim],dz_sq[slicesPerDim];

    // pre-calculate d[xyz]_sq for efficiency
    // 子采样点在各方向的偏移为 -0.4, -0.2, 0, 0.2, 0.4（即每片的中心），
    // 先各自算好平方，三重循环里只做加法和比较。
    for(plint i=0;i<slicesPerDim;i++){
      T const delta = -0.5 + ((T)i+0.5)*sliceWidth;
      T const dx = dx_+delta; dx_sq[i] = dx*dx;
      T const dy = dy_+delta; dy_sq[i] = dy*dy;
      T const dz = dz_+delta; dz_sq[i] = dz*dz;
    }

    // 计数：有多少个子采样点落在球内（布尔值直接当 0/1 相加）
    pluint n(0);
    for(plint i=0;i<slicesPerDim;i++){
      for(plint j=0;j<slicesPerDim;j++){
        for(plint k=0;k<slicesPerDim;k++){
          n += (dx_sq[i] + dy_sq[j] + dz_sq[k] < r_sq);
        }
      }
    }

    // 固含率 = 球内子格数 / 125
    return fraction*((T)n);
  }



  /*
   * setValues：把本颗粒的信息写入格点。
   * 格点处的颗粒速度 = 质心平动速度 v + ω × r，r = (dx,dy,dz) 为格点相对质心的偏移。
   * 叉积按分量展开：
   *   (ω × r)_x = ω_y·dz − ω_z·dy
   *   (ω × r)_y = ω_z·dx − ω_x·dz
   *   (ω × r)_z = ω_x·dy − ω_y·dx
   */
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setValues(IBdynamicsParticleData<T,Descriptor> &p,
                                                  T const sf, T const dx, T const dy, T const dz)
  {
    p.uPart.from_cArray(v);
    if(omega != 0){ // omega 指针可能为空（例如不考虑转动时）
      p.uPart[0] += omega[1]*dz - omega[2]*dy;
      p.uPart[1] += -omega[0]*dz + omega[2]*dx;
      p.uPart[2] += omega[0]*dy - omega[1]*dx;
    }
    p.solidFraction = sf;
    p.partId = id;
  }

  // setToZero：标记该格点为"无颗粒"（partId=0 是关键，SumForceTorque3D 靠它跳过）
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setToZero(IBdynamicsParticleData<T,Descriptor> &p)
  {
    p.uPart[0] = 0;
    p.uPart[1] = 0;
    p.uPart[2] = 0;
    p.solidFraction = 0;
    p.partId = 0;
  }

  /*
   * implementation of SumForceTorque3D
   */

  /* --------------------------------------------- */
  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(typename ParticleData<T>::ParticleDataArrayVector &x_,
                                                   T *force_, T *torque_, LiggghtsCouplingWrapper &wrapper_)
    : x(x_), force(force_), torque(torque_), wrapper(wrapper_)
  {}

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig)
    : BoxProcessingFunctional3D_L<T,Descriptor>(orig), x(orig.x),
      force(orig.force), torque(orig.torque), wrapper(orig.wrapper) {}

  /*
   * process：遍历本地块，把每个被颗粒覆盖格点的水动力累加到对应颗粒。
   * 只在本进程的 bulk 区域（不含幽灵层）上施加，保证每个格点的力只被计一次。
   */
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();

    // "real" domain size is nx-2 etc
    // 这里的 lattice 是本地 BlockLattice3D，其 getNx() 包含两侧各 1 层幽灵层，
    // 所以减 2 才是真实块尺寸；它只在下面的最小镜像判断里用到。
    plint nx = lattice.getNx()-2, ny = lattice.getNy()-2, nz = lattice.getNz()-2;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);

          IBdynamicsParticleData<T,Descriptor>* particleData =
            getParticleDataFromCell<T,Descriptor>(cell);

          if(!particleData) continue;

          // LIGGGHTS indices start at 1
          // partId 是 LIGGGHTS 的全局编号 tag（从 1 起），0 表示此格点无颗粒
          plint const id = particleData->partId;
          if(id < 1) continue; // no particle here

          // 全局编号 → 本进程的本地索引（包含幽灵颗粒），x/force/torque 都按本地索引排列
          plint const ind = wrapper.lmp->atom->map(id);

          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);

          // 力臂 r = 格点位置 − 颗粒中心（格子单位）
	  T dx = xGlobal - x[ind][0];
	  T dy = yGlobal - x[ind][1];
	  T dz = zGlobal - x[ind][2];

          // minimum image convention, needed if
          // (1) PBC are used and
          // (2) both ends of PBC lie on the same processor
          //
          // 最小镜像约定（minimum image convention）：周期边界下，一个颗粒可能
          // "一半在域左端、一半在域右端"，此时格点与颗粒中心的直接差值可能接近
          // 整个域长。若 |dx| 超过半个域长，就把它平移一个域长，取最近的那个像。
          // 注意这里用的是本地块尺寸 nx 而非全局尺寸——只有当周期边界的两端
          // 都落在同一进程（例如单进程或该方向未分块）时这个修正才会触发。
          if(dx > nx/2) dx -= nx;
          else if(dx < -nx/2) dx += nx;
          if(dy > ny/2) dy -= ny;
          else if(dy < -ny/2) dy += ny;
          if(dz > nz/2) dz -= nz;
          else if(dz < -nz/2) dz += nz;

          // 该格点在碰撞步里通过动量交换算出的力（见 ibCompositeDynamics3D.hh）
          T const forceX = particleData->hydrodynamicForce[0];
          T const forceY = particleData->hydrodynamicForce[1];
          T const forceZ = particleData->hydrodynamicForce[2];

          // 力矩 τ = r × F，按分量展开
          T const torqueX = dy*forceZ - dz*forceY;
          T const torqueY = -dx*forceZ + dz*forceX;
          T const torqueZ = dx*forceY - dy*forceX;

          // 累加到该颗粒（多个格点的贡献在这里求和成合力/合力矩）
          addForce(ind,0,forceX);
          addForce(ind,1,forceY);
          addForce(ind,2,forceZ);

          addTorque(ind,0,torqueX);
          addTorque(ind,1,torqueY);
          addTorque(ind,2,torqueZ);
        }
      }
    }
  }


  // force/torque 是扁平数组：第 partId 个颗粒的第 coord 分量放在 3*partId+coord
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addForce(plint const partId, plint const coord, T const value)
  {
    force[3*partId+coord] += value;
  }
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addTorque(plint const partId, plint const coord, T const value)
  {
    torque[3*partId+coord] += value;
  }

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>* SumForceTorque3D<T,Descriptor>::clone() const
  {
    return new SumForceTorque3D<T,Descriptor>(*this);
  }

  // 只读格子，不修改任何格点数据，因此不需要触发幽灵层同步
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
  }


};

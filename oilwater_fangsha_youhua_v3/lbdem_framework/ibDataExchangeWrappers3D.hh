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
 * implementation of data exchange wrapper functions
 */

/*
 * ============================ 文件说明 ============================
 * ibDataExchangeWrappers3D.hh —— setSpheresOnLattice / getForcesFromLattice 的实现，
 * 是 LIGGGHTS 与 Palabos 之间"搬数据"的总调度。
 *
 * 数据流一览：
 *
 *  setSpheresOnLattice（颗粒 → 格子）
 *    LIGGGHTS: atom->x/v/omega/radius（物理单位，本地+幽灵颗粒）
 *      --units.getLb*()--> 格子单位
 *      --new SetSingleSphere3D--> 对每个颗粒，在"包围盒 ∩ 本地块 (∩ 耦合域)"上
 *                                  施加 functional，写入格点 particleData
 *      --AttributeFunctional--> 一次性同步各进程的幽灵层
 *
 *  getForcesFromLattice（格子 → 颗粒）
 *    格点 particleData.hydrodynamicForce
 *      --SumForceTorque3D（只在本地块）--> force[]/torque[]（格子单位，按本地颗粒索引）
 *      --units.getPhysForce/Torque()--> 物理单位
 *      --写入 fix couple/lb/onetoone 的 f/t 数组--> couplingFix->comm_force_torque()
 *        做 MPI 归约，把分散在各进程（含幽灵副本）的贡献加到颗粒所属进程上
 *
 * 为什么要处理"本地 + 幽灵"颗粒（nlocal + nghost）？
 *   LIGGGHTS 按空间分块，每个进程只"拥有"（own）落在自己子域里的颗粒（local），
 *   但靠近边界的颗粒会以"幽灵"（ghost）副本的形式出现在邻居进程。一个颗粒
 *   的球体可能跨越进程边界，它在邻居进程那侧覆盖的格点只能由邻居进程来
 *   映射和计算受力。因此两个方向都必须把幽灵颗粒算进去，最后靠
 *   comm_force_torque 把幽灵副本上累加的力归约回拥有者。
 *
 * 两个前提：
 *   - Palabos 的格子分块与 LIGGGHTS 的域分解一致（latticeDecomposition.cpp 保证），
 *     且每个进程恰好一个块（代码里多处依赖 getBlocks()[0]）；
 *   - LIGGGHTS 输入脚本里定义了 fix couple/lb/onetoone。
 * ================================================================
 */

#ifndef IBDATAEXCHANGEWRAPPERS_HH_LBDEM
#define IBDATAEXCHANGEWRAPPERS_HH_LBDEM

// 多球刚体（multisphere）支持只在 LIGGGHTS-PREMIUM 中有，公开版默认不编译这段
#ifdef LBDEM_USE_MULTISPHERE
#include "fix_multisphere.h"
#endif

namespace plb {
  /*
   * setSpheresOnLattice（主实现）
   * 注意：最后一个参数 couplingDomain 在这里给了默认值 0，所以 5 参版本同时
   * 覆盖了"不限域"和"限域"两种用法；另外两个重载只是转发到这里。
   */
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           std::vector<plint> &excludeType,
                           bool initVelFlag,
                           Box3D const* couplingDomain = 0)
  {
    plint rank = global::mpi().getRank();

    // this relies on the fact that there is exactly one block on each lattice
    // 取本进程唯一的那个块在全局格子坐标中的范围 localBB（bulk 部分，不含幽灵层）
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];

    plint nx=lattice.getNx(), ny=lattice.getNy(), nz=lattice.getNz();
    // 本进程可见的颗粒数 = 自己拥有的 + 幽灵副本（理由见文件说明）
    plint nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;


      // LIGGGHTS-public does not support multisphere (yet?)
      // thus need to explicitly exclude multisphere method
#ifdef LBDEM_USE_MULTISPHERE

    FixMultisphere *fixMultisphere = 0;
    plint nms = wrapper.lmp->modify->n_fixes_style("multisphere");
    if(nms)
      fixMultisphere = static_cast<FixMultisphere*>(wrapper.lmp->modify->find_fix_style("multisphere",0));

#endif /* LBDEM_USE_MULTISPHERE */

    // 逐颗粒处理：每个颗粒构造一个 SetSingleSphere3D 并施加到它影响的区域
    for(plint iS=0;iS<nPart;iS++){
      // 第 1 步：按颗粒类型过滤（excludeType 里列出的类型不参与流固耦合）
      plint type = (plint)wrapper.lmp->atom->type[iS];
      bool excludeFlag(false);
      for(plint iT=0;iT<excludeType.size();iT++){
        //        pcout << iS << " " << type << " " << excludeType[iT] << std::endl;
        if(type == excludeType[iT]){
          excludeFlag = true;
          break;
        }
      }

      if(excludeFlag) continue;

      // 第 2 步：把颗粒的位置、速度、角速度、半径从物理单位换成格子单位
      T x[3],v[3],omega[3];
      T r;
      plint id = wrapper.lmp->atom->tag[iS]; // 全局编号（跨进程唯一，从 1 起）

      SetSingleSphere3D<T,Descriptor> *sss = 0;

      for(plint i=0;i<3;i++){
        x[i] = units.getLbPosition(wrapper.lmp->atom->x[iS][i]);
      }
      r = units.getLbLength(wrapper.lmp->atom->radius[iS]);

      // LIGGGHTS-public does not support multisphere (yet?)
      // thus need to explicitly exclude multisphere method
#ifdef LBDEM_USE_MULTISPHERE

      plint msBody = fixMultisphere ? fixMultisphere->belongs_to(id) : -1;


      // if particle belongs to multisphere, don't use omega
      // 多球刚体：各子球的 v 已经包含了整体转动的贡献，故 omega 置 0 避免重复计入
      if(msBody > -1){
        T x_com[3],v_com[3];
	fixMultisphere->data().xcm(x_com,msBody);
	fixMultisphere->data().vcm(v_com,msBody);
        for(plint i=0;i<3;i++){
          v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
	  omega[i] = 0;//units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
          x_com[i] = units.getLbPosition(x_com[i]);
          v_com[i] = units.getLbVel(v_com[i]);
        }

	sss = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);

      } else {
        for(plint i=0;i<3;i++){
          v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
	  omega[i] = units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
        }

        // use sphere center as center of mass for simple spheres
        sss  = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);
      }

#else
      // 公开版 LIGGGHTS（本案例）：单球颗粒，球心即质心，ω 直接取 atom->omega
      for(plint i=0;i<3;i++){
        v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
        omega[i] = units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
      }
      r = units.getLbLength(wrapper.lmp->atom->radius[iS]);
      sss = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);

#endif /* LBDEM_USE_MULTISPHERE */

      // 第 3 步：确定要施加 functional 的区域 = 颗粒包围盒 ∩ 本进程块
      Box3D sss_box = sss->getBoundingBox();

      // only go over part that lies on local processor
      // to avoid unnecessary communication overhead
      Box3D sss_box_intersect(0,0,0,0,0,0);

      bool boxes_intersect = intersect(sss_box,localBB,sss_box_intersect);

      // 若给定耦合域，再把固含率映射范围限制在耦合域内（全局格子坐标）
      // 本案例用它把出口压力边界附近几层格子排除在外：那里的格点若被标成
      // 固体，会与 Zou-He 压力边界条件互相干扰。
      if(boxes_intersect && couplingDomain){
        Box3D sss_box_coupled(0,0,0,0,0,0);
        if(intersect(sss_box_intersect,*couplingDomain,sss_box_coupled)){
          sss_box_intersect = sss_box_coupled;
        } else {
          boxes_intersect = false;
        }
      }

      // 第 4 步：有交集才真正施加（applyProcessingFunctional 会接管并释放 sss），
      // 没交集说明该颗粒与本进程无关，自己释放掉。
      if(boxes_intersect)
        applyProcessingFunctional(sss,sss_box_intersect,lattice);
      else
        delete sss;
      // applyProcessingFunctional(sss,sss_box,lattice);
    }

    // this one returns modif::staticVariables and forces an update of those along processor
    // boundaries
    // 所有颗粒放完后，施加一个"空"functional 触发一次幽灵层同步：SetSingleSphere3D
    // 为省通信故意声明自己不修改数据，所以 particleData 的跨进程一致性靠这一步保证。
    applyProcessingFunctional(new AttributeFunctional<T,Descriptor>(),lattice.getBoundingBox(),lattice);

  }

  // 4 参重载：不排除类型、不限域 → 转发给主实现
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag)
  {
    std::vector<plint> dummyExcludeType;
    setSpheresOnLattice(lattice,wrapper,units,dummyExcludeType,initVelFlag,0);
  }

  // 带耦合域的重载（本案例新增并实际使用）：不排除类型、限域 → 转发给主实现
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag,
                           Box3D const* couplingDomain)
  {
    std::vector<plint> dummyExcludeType;
    setSpheresOnLattice(lattice,wrapper,units,dummyExcludeType,initVelFlag,couplingDomain);
  }



  /*
   * getForcesFromLattice：格子 → 颗粒
   * 步骤：
   *   1. 准备颗粒位置数组 x_lb（格子单位），SumForceTorque3D 算力臂时要用；
   *   2. 准备并清零 force/torque 累加数组（长度 3*nPart）；
   *   3. 在本地块上施加 SumForceTorque3D，得到每个本地/幽灵颗粒的合力合力矩；
   *   4. 换成物理单位，写进 LIGGGHTS 的 fix couple/lb/onetoone；
   *   5. comm_force_torque() 做 MPI 归约。
   */
  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units)
  {
    // debug stuff
    plint r = global::mpi().getRank();

    // 声明为 static 是为了跨调用复用内存：颗粒数通常不变，避免每步重新分配
    static std::vector<T> force,torque;
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;

    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    plint const n_force = nPart*3;

    if(nPart == 0) return; // no particles - no work

    // 步骤 1：更新 x_lb。若颗粒数比上次多就追加，否则原地覆盖
    //（幽灵颗粒数每步都可能变化，所以不能假设大小固定）
    if(nPart > x_lb.size()){
      for(plint iPart=0;iPart<x_lb.size();iPart++){
        x_lb[iPart][0] = units.getLbPosition(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbPosition(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbPosition(wrapper.lmp->atom->x[iPart][2]);
      }
      for(plint iPart = x_lb.size();iPart < nPart; iPart++)
        x_lb.push_back( Array<T,3>( units.getLbPosition(wrapper.lmp->atom->x[iPart][0]),
                                    units.getLbPosition(wrapper.lmp->atom->x[iPart][1]),
                                    units.getLbPosition(wrapper.lmp->atom->x[iPart][2]) ) );
    } else{
      for(plint iPart=0;iPart<nPart;iPart++){
        x_lb[iPart][0] = units.getLbPosition(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbPosition(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbPosition(wrapper.lmp->atom->x[iPart][2]);
      }
    }

    // 步骤 2：force/torque 清零（同样的"够大就覆盖、不够就追加"逻辑）
    if(n_force > force.size()){
      for(plint i=0;i<force.size();i++){
        force[i] = 0;
        torque[i] = 0;
      }
      for(plint i=force.size();i<n_force;i++){
        force.push_back(0.);
        torque.push_back(0.);
      }
    } else {
      for(plint i=0;i<n_force;i++){
        force[i] = 0;
        torque[i] = 0;
      }
    }

    // 步骤 3：构造求和 functional，并只在本进程的 bulk 块上施加，
    // 保证每个格点的力恰好被计一次（幽灵层格点由拥有它的进程负责）
    SumForceTorque3D<T,Descriptor> *sft = new SumForceTorque3D<T,Descriptor>(x_lb,
                                                                             &force.front(),&torque.front(),
                                                                             wrapper);

    // this relies on the fact that there is exactly one block on each processor
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];
    applyProcessingFunctional(sft,localBB, lattice);


    // 步骤 4：找到 LIGGGHTS 脚本里定义的 fix couple/lb/onetoone，
    // 它在 LIGGGHTS 侧持有 per-atom 的流体力/力矩数组，并在 DEM 积分时把它们
    // 加到颗粒上（"onetoone" 指每个 LBM 颗粒一一对应一个 DEM 颗粒）
    LAMMPS_NS::FixLbCouplingOnetoone
      *couplingFix
      = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
      (wrapper.lmp->modify->find_fix_style("couple/lb/onetoone",0));

    double **f_liggghts = couplingFix->get_force_ptr();
    double **t_liggghts = couplingFix->get_torque_ptr();

    // 先把 LIGGGHTS 侧的力数组清零（上一个耦合步的值不能残留）
    for(plint iPart=0;iPart<nPart;iPart++)
      for(plint j=0;j<3;j++){
        f_liggghts[iPart][j] = 0;
        t_liggghts[iPart][j] = 0;
      }


    // 把格子单位的力/力矩换成物理单位（因子 ρΔx⁴/Δt² 与 ρΔx⁵/Δt²，见 physunits.h）
    // 并按 LIGGGHTS 的本地索引写入。tag→map 这一步在当前排列下一般是恒等映射，
    // 写成这样是为了稳妥。
    for(plint iPart=0;iPart<nPart;iPart++){
      int tag = wrapper.lmp->atom->tag[iPart];
      int liggghts_ind = wrapper.lmp->atom->map(tag);

      for(plint j=0;j<3;j++){
        f_liggghts[liggghts_ind][j] += units.getPhysForce(force[3*iPart+j]);
        t_liggghts[liggghts_ind][j] += units.getPhysTorque(torque[3*iPart+j]);
      }
    }
    // 步骤 5：MPI 归约。同一个颗粒在邻居进程上以幽灵副本形式也累加了力，
    // comm_force_torque 执行"反向通信"（reverse comm）把幽灵副本上的贡献
    // 加回到拥有该颗粒的进程，然后再正向广播，使各进程看到完整的合力。
    couplingFix->comm_force_torque();
  }

}; /* namespace plb */

#endif /* IBDATAEXCHANGEWRAPPERS_HH_LBDEM */

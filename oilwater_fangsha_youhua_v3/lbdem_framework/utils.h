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
 * utils.h —— 一个工具函数：getParticleDataFromCell()。
 *
 * 问题背景：Palabos 中每个格点（Cell）持有一个"动力学对象"（Dynamics），
 * 而且动力学可以层层嵌套（CompositeDynamics 包着 baseDynamics，baseDynamics
 * 还可能再是一个 Composite）。本框架的 IBcompositeDynamics 就是这样一层
 * "外壳"，颗粒信息 particleData 就存在它里面。
 *
 * 所以，当某个 data processor（如 SetSingleSphere3D、SumForceTorque3D、
 * GetScalarQuantityFromDynamicsFunctional）想读写某格点的颗粒信息时，
 * 就需要一个办法"剥开动力学的洋葱皮，找到 IBcompositeDynamics 那一层"。
 * getParticleDataFromCell() 正是做这件事：找到就返回 particleData 的指针，
 * 找不到（该格点不是 IB 动力学，例如边界节点或 NoDynamics）就返回 0。
 *
 * 被谁使用：ibProcessors3D.hh、ibDataWritingFunctionals3D.hh。
 * ================================================================
 */


#ifndef UTILS_H_LBDEM
#define UTILS_H_LBDEM

#include "ibDynamicsParticleData.h"

namespace plb {
  template<typename T, template<typename U> class Descriptor>
  IBdynamicsParticleData<T,Descriptor>* getParticleDataFromCell(Cell<T,Descriptor> &cell);


}; /* namespace plb */

// --------------------------
// implementaton starts here

namespace plb {
  template<typename T, template<typename U> class Descriptor>
  IBdynamicsParticleData<T,Descriptor>* getParticleDataFromCell(Cell<T,Descriptor> &cell)
  {
    // 第 1 步：拿到 IBcompositeDynamics 在 Palabos 动力学注册表中的数字 id。
    // 用 static 构造一个"样板"对象只为了读它的 id，整个程序只执行一次，
    // 之后每次调用都直接比较整数，开销极小。
    static IBcompositeDynamics<T,Descriptor> const cdyn(new NoDynamics<T,Descriptor>(),false);
    static plint const cdynId = cdyn.getId();

    Dynamics<T,Descriptor> *dyn = &( cell.getDynamics() );

    // 第 2 步：最外层就是 IB 动力学的常见情况，直接返回。
    if(dyn->getId() == cdynId){
      return &( (static_cast< IBcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }

    // 第 3 步：否则逐层向内剥，只要当前层是 Composite 就取它的 baseDynamics 再比较。
    // 这样即使 IB 动力学外面又套了别的 Composite（例如 Smagorinsky 包装），也能找到。
    while(dyn->isComposite()){
      dyn = &(static_cast<CompositeDynamics<T,Descriptor>* >(dyn))->getBaseDynamics();
      if(dyn->getId() == cdynId)
        return &( (static_cast< IBcompositeDynamics<T,Descriptor>* >(dyn))->particleData );
    }

    // 第 4 步：剥到底都没有 IB 动力学——该格点不参与颗粒耦合。
    return 0;
  }

}; /* namespace plb */

#endif /* UTILS_H_LBDEM */

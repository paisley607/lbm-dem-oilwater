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
 * composite dynamics that implement the method of Noble and Torczinsky
 */

/*
 * ============================ 文件说明 ============================
 * ibCompositeDynamics3D.h —— PSM/IBM 碰撞算子的类声明（流体侧的核心）。
 *
 * IBcompositeDynamics 是 Palabos 的 CompositeDynamics（复合动力学）的子类：
 * 它"包住"一个普通的流体碰撞算子（baseDynamics，例如 BGK 或本案例使用的
 * Smagorinsky 正则化动力学），在普通碰撞之后，再按本格点的固含率
 * （solid fraction）ε 叠加 Noble–Torczynski 的固体修正项，并顺便算出
 * 该格点对颗粒的水动力（hydrodynamic force）。
 *
 * 关键成员：
 *   particleData —— 每个格点一份的颗粒信息（ε、颗粒速度、颗粒 id、受力），
 *                   见 ibDynamicsParticleData.h；
 *   collide()    —— 碰撞步，PSM 公式在此实现（见 .hh）；
 *   prepareCollision() —— 碰撞前备份碰撞前分布 fPre，PSM 公式需要它；
 *   defineVelocity()   —— 用给定速度把格点初始化到平衡态（颗粒内部初始化用）；
 *   serialize/unserialize —— Palabos 跨进程传递动力学对象时打包/解包 particleData。
 *
 * 与其它文件的关系：
 *   主程序通过 new IBcompositeDynamics(new 某流体动力学) 给整个格子赋动力学；
 *   ibProcessors3D.hh 写 particleData（映射固含率）、读 particleData（汇总受力）；
 *   utils.h 提供从 Cell 找到本类 particleData 的工具。
 * ================================================================
 */

#ifndef IB_COMPOSITE_DYNAMICS_H_LBDEM
#define IB_COMPOSITE_DYNAMICS_H_LBDEM

#include "ibDynamicsParticleData.h"

namespace plb {

  template<typename T, template<typename U> class Descriptor>
  class IBcompositeDynamics : public CompositeDynamics<T,Descriptor>
  {
  public:

    // baseDynamics_：被包装的普通流体碰撞算子（所有权转移给本对象）。
    // automaticPrepareCollision_：是否由 Palabos 在 collide 前自动调 prepareCollision；
    // 本类的 collide() 内部会自己调一次，所以主程序一般用默认值即可。
    IBcompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                        bool automaticPrepareCollision_ = true);
    IBcompositeDynamics(const IBcompositeDynamics &orig);
    // 反序列化构造：Palabos 在进程间迁移/重建格点时用它从字节流恢复对象
    IBcompositeDynamics(HierarchicUnserializer &unserializer);
    ~IBcompositeDynamics();

    virtual IBcompositeDynamics<T,Descriptor>* clone() const;

    // 返回本类在 Palabos 动力学注册表中的 id（utils.h 用它来识别 IB 格点）
    virtual int getId() const;

    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);

    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    virtual void collide(Cell<T,Descriptor>& cell,  BlockStatistics& statistics);

    virtual void defineVelocity(Cell<T,Descriptor>& cell,
                                Array<T,Descriptor<T>::d> const& u);

    // 公开成员：本格点的颗粒信息，由 data processor 直接读写
    IBdynamicsParticleData<T,Descriptor> particleData;
  private:
    static int id;

    // performance tweaks
    // // have static variables for some temporary quantities
    // // to save memory
    //
    // 以下三个临时数组声明为 static（全类共享一份），而不是每个格点各存一份：
    // 对几百万个格点来说能省下大量内存。它们只在单个 collide() 调用内有效，
    // 调用结束即失去意义（因此本类不是线程安全的，Palabos 的 MPI 并行不受影响）。
    static Array<T,Descriptor<T>::q> fEqSolid; // 以颗粒速度 uPart 计算的平衡分布 f_eq(ρ, u_s)
    static Array<T,Descriptor<T>::q> fEq;      // 以流体速度 u 计算的平衡分布 f_eq(ρ, u)
    static Array<T,Descriptor<T>::q> fPre;     // 碰撞前的分布函数备份 f_pre
  };

}; /* namespace plb */

#include "ibCompositeDynamics3D.hh"

#endif /* IB_COMPOSITE_DYNAMICS3D_H_LBDEM */

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
 * ibDynamicsParticleData.h —— "每个格点上挂一份的颗粒信息"的数据结构。
 *
 * 在 PSM/IBM 耦合中，流体格点（lattice node）需要知道"我被哪个颗粒盖住了、
 * 盖住了多少、颗粒在我这里的速度是多少"，碰撞步才能据此做加权修正，
 * 并把动量交换得到的力记录下来。本文件定义的 IBdynamicsParticleData 就是
 * 这份"随格点走"的小结构体，它作为成员变量嵌在 IBcompositeDynamics 里
 * （见 ibCompositeDynamics3D.h），因此每个使用 IB 动力学的格点都各自持有一份。
 *
 * 数据流：
 *   写入：ibProcessors3D.hh 的 SetSingleSphere3D 在每个耦合步把 LIGGGHTS 的
 *         颗粒位置/速度映射成 solidFraction、uPart、partId；
 *   更新：ibCompositeDynamics3D.hh 的 collide() 在碰撞时算出 hydrodynamicForce；
 *   读出：SumForceTorque3D 把 hydrodynamicForce 按 partId 汇总回颗粒；
 *         ibDataWritingFunctionals3D 把这些量导出成场数据用于后处理。
 *
 * 配套的 .hh 文件为空，仅保留文件结构的对称性。
 * ================================================================
 */

#ifndef IB_DYNAMICS_PARTICLE_DATA_H_LBDEM
#define IB_DYNAMICS_PARTICLE_DATA_H_LBDEM

namespace plb {

  // T：浮点类型（一般为 double）；Descriptor：格子模型（如 D3Q19），
  // Descriptor<T>::d 为空间维数（3）。
  template<typename T, template<typename U> class Descriptor>
  struct IBdynamicsParticleData {
  public:
    // 默认构造：无颗粒覆盖——id 为 0（LIGGGHTS 的颗粒编号从 1 开始，0 表示"没有颗粒"），
    // 固含率为 0，颗粒速度和水动力为零向量。
    IBdynamicsParticleData()
      : partId(0), solidFraction(0.),
        uPart(0.,0.,0.), hydrodynamicForce(0.,0.,0.)
    {}
    plint partId;          // 覆盖本格点的颗粒在 LIGGGHTS 中的全局编号（tag），0 = 无颗粒
    T solidFraction;       // 固含率（solid fraction）ε：本格点单位立方体被颗粒占据的体积比例，0~1
    Array<T,Descriptor<T>::d> uPart;             // 颗粒在本格点处的速度（格子单位），
                                                 // = 质心平动速度 + ω × r，r 为格点相对质心的位置
    Array<T,Descriptor<T>::d> hydrodynamicForce; // 本格点在最近一次碰撞中通过动量交换
                                                 // 传递给颗粒的力（格子单位），每步重新清零累加
  };

}; /* namespace plb */


#endif /* IB_DYNAMICS_PARTICLE_DATA_H_LBDEM */

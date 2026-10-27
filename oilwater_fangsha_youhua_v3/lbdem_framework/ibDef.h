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
 * ibDef.h —— LBM-DEM 耦合框架的全局宏定义（只有 3 个宏，但影响全局行为）。
 *
 * 本框架采用 Noble–Torczynski 的部分饱和法（Partially Saturated Method, PSM），
 * 代码里历史上叫 "IB"（immersed boundary），两者在此指同一套方法。
 * 核心思想：每个流体格点（lattice node）记录一个"固含率"（solid fraction，
 * 记作 ε 或 sf，取值 0~1），表示该格点所在的单位立方体被颗粒占据的体积比例。
 * 碰撞步按 ε 在"普通流体碰撞"和"固体反弹"之间加权，颗粒受力由动量交换累加。
 *
 * 三个宏：
 *   LBDEM_USE_WEIGHING —— 是否使用带松弛时间 τ 的 B 加权函数（见 ibCompositeDynamics3D.hh）。
 *   SOLFRAC_MIN / SOLFRAC_MAX —— 判断"格点是否被颗粒覆盖/完全覆盖"的阈值。
 *
 * 被谁使用：ibCompositeDynamics3D.hh（碰撞）、ibProcessors3D.hh（固含率映射）。
 * ================================================================
 */

#ifndef IB_DEF_H
#define IB_DEF_H

/*
 * LBDEM_USE_WEIGHING：开启后，碰撞步中的固体权重 B 不再直接等于固含率 ε，
 * 而是 B = ε(τ-0.5) / ((1-ε) + (τ-0.5))，τ 为 LBM 松弛时间。
 * 这是 Noble & Torczynski (1998) 原文推荐的形式，能减弱 B 对 τ 的依赖，
 * 使不同粘度下的受力更一致。
 *
 * 注意：这个宏是在 "#ifdef" 里被检查的，而 ibCompositeDynamics3D.hh 只有在
 * 包含了 ibDef.h 之后才会看到它。主程序（oilwaterDrive.cpp）在包含任何框架
 * 头文件之前又自己 #define 了一次，目的是无论头文件包含顺序如何、无论编译
 * 命令是否传 -DLBDEM_USE_WEIGHING，都保证主程序里实例化的模板和框架内部
 * 使用同一套加权公式；若两处不一致，同一个模板可能被编成两种行为（ODR 违规）。
 */
#define LBDEM_USE_WEIGHING

/*
 * SOLFRAC_MIN：固含率低于此值的格点视为"纯流体"，跳过 PSM 修正，直接用普通碰撞。
 * SOLFRAC_MAX：固含率高于此值的格点视为"完全在颗粒内部"，跳过流体碰撞，
 *              直接采用完全反弹（bounce-back）形式的非平衡反弹。
 * 介于两者之间的格点才做 ε 加权的混合碰撞。
 * 取 0.001/0.999 而不是 0/1 是为了吸收 5x5x5 子采样和浮点运算带来的舍入误差。
 */
#define SOLFRAC_MIN 0.001
#define SOLFRAC_MAX 0.999

#endif // IB_DEF_H

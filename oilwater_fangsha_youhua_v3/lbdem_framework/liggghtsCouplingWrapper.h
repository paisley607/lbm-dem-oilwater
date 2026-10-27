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
 * liggghtsCouplingWrapper.h —— LiggghtsCouplingWrapper 类的声明：
 * 把 LIGGGHTS（DEM 求解器）当成一个"库"嵌进 Palabos 主程序里的薄封装。
 *
 * LIGGGHTS 本身是个独立程序，靠读脚本运行；"库模式"（library mode）允许
 * 其它 C++ 程序 new 一个 LAMMPS 对象，然后逐条喂它命令。本类就是这么做的：
 *   - 构造函数：在给定的 MPI 通信子上创建 LAMMPS 实例；
 *   - execFile：执行一整个 LIGGGHTS 输入脚本（如 in.lbdem，定义颗粒、接触模型、
 *               fix couple/lb/onetoone 等）；
 *   - execCommand：执行单条命令；
 *   - setVariable：在脚本里定义/覆盖一个 variable（把 C++ 侧参数传给脚本）；
 *   - run / runUpto / runWithoutSetup：推进 DEM 若干步。
 *
 * 公开成员 lmp 是 LAMMPS 对象本身，耦合代码直接通过它访问颗粒数组
 * （lmp->atom->x, v, omega, radius, tag, nlocal, nghost）和 fix 列表（lmp->modify）。
 *
 * 实现见 liggghtsCouplingWrapper.cpp（非模板，需单独编译链接）。
 * 调用者：主程序、ibDataExchangeWrappers3D.hh、ibProcessors3D.hh、latticeDecomposition。
 * ================================================================
 */

#ifndef LIGGGHTSCOUPLINGWRAPPER_H
#define LIGGGHTSCOUPLINGWRAPPER_H

#include "palabos3D.h"
#include "palabos3D.hh"

// necessary LAMMPS/LIGGGHTS includes

#include "lammps.h"              // LAMMPS_NS::LAMMPS 类
#include "input.h"               // lmp->input->file() / one()
#include "library.h"             // C 接口：lammps_get_natoms 等
#include "library_cfd_coupling.h" // LIGGGHTS 为 CFD 耦合提供的额外库接口

class LiggghtsCouplingWrapper {
public:
  // argv：主程序的命令行（只用 argv[0] 作为程序名传给 LAMMPS）；
  // communicator：LIGGGHTS 使用的 MPI 通信子，与 Palabos 共用同一个即可
  LiggghtsCouplingWrapper(char **argv, MPI_Comm communicator);
  void execFile(char* const fname);                  // 执行整个脚本文件
  void execCommand(std::stringstream const &cmd);    // 执行单条命令（stringstream 版）
  void execCommand(char* const cmd);                 // 执行单条命令（C 字符串版）
  void run(plb::plint nSteps);                       // "run N"：再推进 N 步
  void runUpto(plb::plint nSteps);                   // "run N upto"：推进到总步数 N
  // 复用已有 DEM setup（拓扑未变时）：等价于 "run nSteps pre no"，
  // 跳过本运行前的邻居/setup 重建。
  void runWithoutSetup(plb::plint nSteps);
  int getNumParticles();                             // 全局颗粒总数（跨进程求和）
  void setVariable(char const *name, double value);  // "variable name equal value"
  void setVariable(char const *name, std::string &value); // "variable name string value"

  //private:
  // 故意不设为 private：耦合代码需要直接读写 LIGGGHTS 内部数据
  LAMMPS_NS::LAMMPS *lmp;
};

#endif /* LIGGGHTSCOUPLINGWRAPPER_H */

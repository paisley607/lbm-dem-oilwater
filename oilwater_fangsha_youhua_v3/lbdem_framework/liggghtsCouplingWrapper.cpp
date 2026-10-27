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
 * liggghtsCouplingWrapper.cpp —— LiggghtsCouplingWrapper 的实现。
 *
 * 所有方法本质上都是"拼一条 LIGGGHTS 命令字符串，交给 lmp->input->one() 执行"，
 * 等价于在 LIGGGHTS 脚本里手写这条命令。几个要点：
 *
 *   execFile      —— 相当于命令行 "liggghts < in.lbdem"，一次性跑完脚本；
 *                    脚本里通常只建模型、不写 run，把时间推进交给主程序控制。
 *   setVariable   —— 把 C++ 里算好的参数（时间步长、流体密度……）注入脚本变量，
 *                    脚本中用 ${name} 引用；必须在用到该变量的命令之前调用。
 *   run(N)        —— 推进 N 个 DEM 步；每次 run 前 LIGGGHTS 会做一次 setup
 *                    （重建邻居表、初始化 fix 等）；
 *   runUpto(N)    —— 推进到全局步数 N（而不是"再走 N 步"），便于与 LBM 步数对齐；
 *   runWithoutSetup(N) —— 本案例新增："run N pre no"，跳过 setup。在耦合循环里
 *                    每个耦合步都要 run 一次，setup 的开销相对很大；颗粒数、
 *                    fix 列表等拓扑不变时可以安全跳过。
 *
 * 典型耦合循环：setSpheresOnLattice → collideAndStream×k → getForcesFromLattice
 *               → wrapper.run(m) 或 runWithoutSetup(m) → 重复。
 * ================================================================
 */

#include "liggghtsCouplingWrapper.h"
#include "mpi.h"

LiggghtsCouplingWrapper::LiggghtsCouplingWrapper(char **argv, MPI_Comm communicator)
  : lmp(0)
{
  // todo: get LAMMPS to recognize command line options
  // 只把程序名 argv[0] 传给 LAMMPS，不转发其它命令行参数
  //（否则 Palabos 的参数会被 LAMMPS 当成自己的选项而报错）
  int argc_lmp = 1;
  char **argv_lmp = 0;
  argv_lmp = new char*[1];
  argv_lmp[0] = argv[0];

  // 在指定通信子上创建 LAMMPS/LIGGGHTS 实例；此后所有进程都持有一个 lmp
  lmp = new LAMMPS_NS::LAMMPS(argc_lmp,argv_lmp,communicator);

  //    delete[] argv_lmp[0];
  // argv_lmp[0] 指向主程序的 argv[0]，不归本类所有，只释放指针数组本身
  delete[] argv_lmp;
}
// 执行整个脚本文件（阻塞直到脚本读完）
void LiggghtsCouplingWrapper::execFile(char* const fname)
{
  lmp->input->file(fname);
}
// 执行单条命令：input->one() 会像读脚本一样解析并执行这一行
void LiggghtsCouplingWrapper::execCommand(std::stringstream const &cmd)
{
  lmp->input->one(cmd.str().c_str());
}
void LiggghtsCouplingWrapper::execCommand(char* const cmd)
{
  lmp->input->one(cmd);
}
// 返回全局颗粒总数（LAMMPS 库接口已做跨进程求和）
int LiggghtsCouplingWrapper::getNumParticles()
{
  return lammps_get_natoms(lmp);
}
// 定义数值型变量：等价于脚本中的 "variable name equal value"
// pcout 只在 rank 0 打印，方便在日志里核对传入的参数
void LiggghtsCouplingWrapper::setVariable(char const *name, double value)
{
  std::stringstream cmd;
  cmd << "variable " << name << " equal " << value;
  plb::pcout << cmd.str() << std::endl;
  execCommand(cmd);
}
// 定义字符串变量：等价于 "variable name string value"（例如传文件名、目录）
void LiggghtsCouplingWrapper::setVariable(char const *name, std::string &value)
{
  std::stringstream cmd;
  cmd << "variable " << name << " string " << value;
  plb::pcout << cmd.str() << std::endl;
  execCommand(cmd);
}
// "run N"：从当前步再推进 N 步，每次调用前 LIGGGHTS 都会执行 setup
void LiggghtsCouplingWrapper::run(plb::plint nSteps)
{
  std::stringstream cmd;
  cmd << "run " << nSteps;
  execCommand(cmd);
}
// "run N upto"：推进到全局时间步 N；若当前已 >= N 则不动
void LiggghtsCouplingWrapper::runUpto(plb::plint nSteps)
{
  std::stringstream cmd;
  cmd << "run " << nSteps << " upto";
  execCommand(cmd);
}
// "run N pre no"：跳过 setup 直接推进；仅当两次 run 之间没有改变颗粒数、
// fix/compute 列表、邻居设置等"拓扑"时才安全（主程序用 demNeedsSetup 标志控制）
void LiggghtsCouplingWrapper::runWithoutSetup(plb::plint nSteps)
{
  // pre no：跳过本次 run 前的 setup（邻居表等），复用上一 run 的状态。
  std::stringstream cmd;
  cmd << "run " << nSteps << " pre no";
  execCommand(cmd);
}





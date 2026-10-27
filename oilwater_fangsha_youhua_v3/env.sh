#!/bin/bash
# ============================================================================
# env.sh —— LBM-DEM 运行环境配置【示例】（北京超算 BSCC-A 分区，账号 sc52079）
#
# 【用途】在登录节点或作业脚本里先执行
#     source /public3/home/sc52079/LiuMingJie/env.sh
#   把编译器、MPI、依赖库的路径加载好，之后再 make 或 bash run_production.sh。
#
# 【换集群要改哪几行】
#   第 ① 组：module 系统的初始化脚本和模块名（每个集群的模块名都不同，用 `module avail` 查）；
#             没有 module 系统的机器可以整段删掉，改成手动 export PATH/CC 等。
#   第 ② 组：LBMDEM_BASE 改成你放 Palabos-PFM、LIGGGHTS-INL 的上级目录；
#             VTK 不是必需的（只有可视化后处理用到），没有就从 LD_LIBRARY_PATH 里删掉。
#   第 ③ 组一般不用动。
# 注意：run_production.sh 对 PALABOS_ROOT/LIGGGHTS_ROOT 另有默认值，env.sh 里 export 的值优先。
# ============================================================================

# ① 编译器 / MPI / 依赖库模块（集群相关）
source /public3/soft/modules/module.sh          # 初始化 module 命令
module purge                                    # 清掉已加载模块，避免版本冲突
module load mpi/openmpi/4.1.1-gcc12.2           # OpenMPI 4.1.1 + GCC 12.2（提供 mpicxx、mpirun）
module load mesa/13.0.6-thc-public3             # Mesa 软件 OpenGL，VTK 离屏渲染用
module load boost/174-gcc-cjj-public3           # Boost 1.74，Palabos 部分头文件依赖
module load miniforge/24.11                     # Python 环境（scons、numpy、后处理脚本）

# ② 代码与库的安装位置（路径相关）
export LBMDEM_BASE=/public3/home/sc52079/LiuMingJie
export PALABOS_ROOT=${LBMDEM_BASE}/Palabos-PFM                    # Palabos（流体 LBM 库）
export LIGGGHTS_ROOT=${LBMDEM_BASE}/LIGGGHTS-INL                  # LIGGGHTS-INL（DEM 库，含平行键）
export LBDEM_ROOT=${LBMDEM_BASE}/LBDEMcoupling-public-master      # 原始耦合框架（本项目改用自带 lbdem_framework，此变量仅供其它案例使用）
export LD_LIBRARY_PATH=${LIGGGHTS_ROOT}/src:${LBMDEM_BASE}/vtk-9.1-install/lib64:${LD_LIBRARY_PATH:-}  # 运行时动态库搜索路径
export PATH=${HOME}/.local/bin:${PATH}                            # 用户级 pip 安装的命令（如 scons）

# ③ 运行时开关
export OMP_NUM_THREADS=1              # 关闭 OpenMP 多线程，只用 MPI 并行，避免超额占用核
export MESA_GL_VERSION_OVERRIDE=3.3   # 让软件 OpenGL 报告 3.3 版本，VTK 才能正常离屏渲染

echo 'LBM-DEM 环境已加载: GCC 12.2 + OpenMPI 4.1.1 + VTK 9.1'

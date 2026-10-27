#!/bin/bash
# =====================================================================================
# submit.sh —— 北京超算 amd_256 分区作业脚本（2026-10-26 修订）
#   提交：sbatch submit.sh        查看：squeue -u $USER        日志：slurm-<作业号>.out
#   取消：scancel <作业号>
#
# 【10-26 改了什么】
#   1. 加了 -t（最长运行时间，单位分钟）。课题组模板里都有这一行；不写时系统按分区默认值（可能很长）
#      计算，容易超过账号的时间限额被拒。冒烟测试估计 2–3 小时，这里给 360 分钟（6 小时）留余量。
#   2. cd 改成 "$SLURM_SUBMIT_DIR"（= 你执行 sbatch 时所在的目录）。上一版写死了旧文件夹
#      .../oilwater_fangsha_youhua，会跑到旧代码去；现在在哪个文件夹里 sbatch 就跑哪个。
#   3. 输出文件名带作业号：slurm-<作业号>.out。
#
# 【核数】amd_256 一个节点 = 64 核。-N（节点数）和 -n（总核数）要配套：-N 1 -n 64，-N 2 -n 128，……
#         改 -n 时，最后一行的 NP 也要改成同一个数。
# =====================================================================================
#SBATCH -J oilwater_smoke
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 600
#SBATCH -o slurm-%j.out

NOW_TIME=`date '+%Y-%m-%d_%H:%M:%S'`
echo "Beginning time is ${NOW_TIME}"

source /public3/home/sc52079/LiuMingJie/env.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
echo "运行目录: $(pwd)"

# ---- ① 冒烟测试（默认）：N=20、只注 0.01 PV，约 2–3 小时 ----
#export CONFIG_FILE=config/test_1pv.env
#python3 plan_run.py --np 64 || echo "[WARN] plan_run.py 认为这组参数通不过预检，程序启动后会报 [FATAL]"
#NP=64 bash run_production.sh

# ---- ② 正式算例：冒烟测试通过、并用实测速度跑过 plan_run.py 之后，
#         把上面三行注释掉，换成下面三行（CONFIG_FILE 换成你要的剖面），并把 -t 改大 ----
 export CONFIG_FILE=config/experiments/run_1pv_m2_n20.env
 python3 plan_run.py --np 64 --mlups-per-core 0.5
 NP=64 bash run_production.sh

NOW_TIME=`date '+%Y-%m-%d_%H:%M:%S'`
echo "End time is ${NOW_TIME}"

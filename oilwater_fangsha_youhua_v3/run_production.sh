#!/bin/bash
# ============================================================================
# run_production.sh —— oilwater_fangsha_youhua 的一键可复现运行脚本
#   物理模型：油-水两相（Shan-Chen）+ 平行胶结细粒 + 出口筛网防砂。
#
# 【整体流程】（下面每一段都有对应的注释标题）
#   ① 读参数   : config/defaults.env 为底，CONFIG_FILE 指定的实验文件覆盖，环境变量优先级最高
#   ② 定运行名 : EXPERIMENT_NAME → 输出目录 OUTDIR（已存在则拒绝，防止覆盖旧结果）
#   ③ 存快照   : 把本次真正用到的关键参数写进 run_metadata.txt
#   ④ 校验几何 : 用 Python 检查颗粒包尺寸、孔隙度、筛网位置是否自洽
#   ⑤ 编译     : make（每次都重编，保证最新头文件进入可执行文件）
#   ⑥ 运行     : mpirun -np NP ./oilwaterDrive 参数...，标准输出存到 run_stdout.log
#   ⑦ 汇总     : 从日志里抓错误/完成/PV 行写 validation_summary.txt，并返回程序退出码
#
# 【用法】
#   bash run_production.sh                                     # 默认参数生产运行
#   CONFIG_FILE=config/experiments/production_v006.env bash run_production.sh   # 指定实验文件
#   RUN_MODE=calibration bash run_production.sh                # 标定初始岩心渗透率（冻结 DEM、全水、低速）
#
# 【常用可覆盖的环境变量】
#   PALABOS_ROOT, LIGGGHTS_ROOT, NP, EXPERIMENT_NAME, OUTDIR, TOTAL_PV,
#   INITIAL_CONTINUOUS_PHASE, INJECT_PHASE, SHAN_CHEN_G, CONFIG_FILE
#   例：NP=8 TOTAL_PV=5 bash run_production.sh
# ============================================================================

# -u：引用未定义变量即报错；-o pipefail：管道中任一命令失败则整条失败
set -uo pipefail

# 案例目录 = 本脚本所在目录（自动定位，不写死路径，整个目录可随意搬动）
CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CASE_DIR" || exit 1

# 依赖库路径：默认指向 HPC 上的安装位置，换机器用环境变量覆盖即可
export PALABOS_ROOT="${PALABOS_ROOT:-/public3/home/sc52079/LiuMingJie/Palabos-PFM}"
export LIGGGHTS_ROOT="${LIGGGHTS_ROOT:-/public3/home/sc52079/LiuMingJie/LIGGGHTS-INL}"
export LD_LIBRARY_PATH="${LIGGGHTS_ROOT}/src:${LD_LIBRARY_PATH:-}"   # 运行时要能找到 libliggghts.so

# ====================== ① 统一参数入口 ======================
# defaults.env 列出全部合法参数；CONFIG_FILE 可指向实验文件（支持 BASE_CONFIG 继承）。
# 这两个路径同时 export 给 C++，主程序用同一套规则再读一遍，保证 bash 与 C++ 看到的参数一致。
PARAMETER_DEFAULTS_FILE="$CASE_DIR/config/defaults.env"
CONFIG_FILE="${CONFIG_FILE:-$PARAMETER_DEFAULTS_FILE}"
export PARAMETER_DEFAULTS_FILE CONFIG_FILE
# shellcheck source=config/load_parameters.sh
source "$CASE_DIR/config/load_parameters.sh" || exit 64
load_water_fangsha_parameters "$PARAMETER_DEFAULTS_FILE" "$CONFIG_FILE" || exit 64

# ====================== ② 运行名 / 输出目录 ======================
# 没给 EXPERIMENT_NAME 就按时间戳+进程号生成一个；名字只允许字母数字和 _.-，避免路径问题
if [[ -z "${EXPERIMENT_NAME:-}" ]]; then
    EXPERIMENT_NAME="oilwater_fangsha_$(date +%Y%m%d_%H%M%S)_${BASHPID}"
fi
if [[ ! "$EXPERIMENT_NAME" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$ ]]; then
    echo "[FATAL] EXPERIMENT_NAME 非法: '$EXPERIMENT_NAME'" >&2
    exit 64
fi
RUN_NAME="$EXPERIMENT_NAME"
OUTPUT_PREFIX="${EXPERIMENT_NAME}__"                    # 输出文件名前缀，便于区分同目录下不同实验
if [[ -z "${OUTDIR:-}" ]]; then OUTDIR="${OUTPUT_ROOT}/${EXPERIMENT_NAME}"; fi

# 标定模式：冻结 DEM、全水饱和、低速注入，只为测初始岩心渗透率。
# 这里用 CALIBRATION_* 参数覆盖对应的生产参数。
CALIBRATION_MODE=0
if [ "${RUN_MODE:-production}" = "calibration" ]; then
    CALIBRATION_MODE=1
    N="$CALIBRATION_N"
    UMAX="$CALIBRATION_UMAX"
    V_INF="$CALIBRATION_V_INF"
    MAX_T="$CALIBRATION_MAX_T"
    RAMP_TIME="$CALIBRATION_RAMP_TIME"
    NU_SCALE=1
    OPEN_OUTLET_PV=0
fi

# 输出目录保护：不能为空或根目录，不能已存在（防止覆盖上次结果）
OUTDIR="${OUTDIR%/}"
if [ -z "$OUTDIR" ] || [ "$OUTDIR" = "/" ]; then
    echo "[FATAL] OUTDIR 为空或为根目录" >&2; exit 64
fi
if [ -e "$OUTDIR" ]; then
    echo "[FATAL] 输出目录已存在，拒绝覆盖: $OUTDIR" >&2; exit 65
fi
if ! mkdir -p "$OUTDIR/tmp" "$OUTDIR/post"; then       # tmp/ 放 LBM 中间文件，post/ 放后处理结果
    echo "[FATAL] 无法创建输出目录: $OUTDIR" >&2; exit 73
fi

# ====================== ③ 保存本次实际使用的参数快照 ======================
# 写成 key=value 文本，后处理脚本和人都能读；dx_um 由参考直径/N 现算
{
    echo "case=$CASE_DIR"
    echo "experiment_name=$EXPERIMENT_NAME"
    echo "captured_at=$(date --iso-8601=seconds)"
    echo "run_mode=${RUN_MODE:-production}"
    echo "calibration_mode=$CALIBRATION_MODE"
    echo "np=$NP"
    echo "N=$N"
    echo "dx_um=$(awk -v d="$REFERENCE_DIAMETER_UM" -v n="$N" 'BEGIN { printf "%.12g", d/n }')"
    echo "total_PV=$TOTAL_PV"
    echo "pack_file=$PACK_FILE"
    echo "rho_o=$RHO_O mu_o=$MU_O shanChenG=$SHAN_CHEN_G"
    echo "initial_continuous_phase=$INITIAL_CONTINUOUS_PHASE inject_phase=$INJECT_PHASE"
    echo "initial_front_position_m=$INITIAL_FRONT_POSITION_M"
    # 本版新增的两相 / 驱动协议参数，一并记录便于事后核对
    echo "tau_w_target=$TAU_W_TARGET"
    echo "wall_contact_angle_deg=$WALL_CONTACT_ANGLE_DEG_TYPE1/$WALL_CONTACT_ANGLE_DEG_TYPE2/$WALL_CONTACT_ANGLE_DEG_TYPE3 container=$CONTAINER_WALL_CONTACT_ANGLE_DEG gain=$WETTING_GAIN"
    echo "inlet_control=$INLET_CONTROL"
    echo "drive_protocol=$DRIVE_PROTOCOL"
} > "$OUTDIR/run_metadata.txt"

echo "======================================================"
echo " oilwater_fangsha_youhua mode=${RUN_MODE:-production} N=$N PV=$TOTAL_PV np=$NP"
echo " rho_o=$RHO_O mu_o=$MU_O G=$SHAN_CHEN_G inject=$INJECT_PHASE"
echo " 输出目录: $OUTDIR"
echo "======================================================"

# ====================== ④ 初始几何校验 ======================
# 几何与流体相无关，复用单相版本的校验脚本。参数单位按脚本要求换成 μm / mm（awk 现算）。
# 校验内容：盒子尺寸、岩心区间、筛网位置与丝径、细粒最小直径、孔隙度范围；结果写 pack_geometry.csv。
if ! python3 validate_water_fangsha_pack.py "$PACK_FILE" \
    --spacing-um "$(awk -v d="$REFERENCE_DIAMETER_UM" -v n="$N" 'BEGIN { print d/n }')" \
    --lx-mm "$(awk -v x="$DOMAIN_LX_M" 'BEGIN { print x*1000 }')" \
    --ly-mm "$(awk -v x="$DOMAIN_LY_M" 'BEGIN { print x*1000 }')" \
    --lz-mm "$(awk -v x="$DOMAIN_LZ_M" 'BEGIN { print x*1000 }')" \
    --core-x0-mm "$(awk -v x="$CORE_X0_M" 'BEGIN { print x*1000 }')" \
    --core-x1-mm "$(awk -v x="$CORE_X1_M" 'BEGIN { print x*1000 }')" \
    --screen-x-mm "$(awk -v x="$SCREEN_X_M" 'BEGIN { print x*1000 }')" \
    --screen-wire-diameter-um "$SCREEN_WIRE_DIAMETER_UM" \
    --fine-diameter-um "$MIN_WEAK_DIAMETER_UM" \
    --min-porosity "$PHI_MIN" --max-porosity "$PHI_MAX" \
    --output "$OUTDIR/pack_geometry.csv"; then
    echo "[FATAL] 初始颗粒包几何校验失败，停止。" >&2
    echo "1" > "$OUTDIR/exit_code.txt"; exit 1
fi

# ====================== ⑤ 编译 ======================
# 每次都重新 make（SCons 只重编有改动的文件，所以不慢），编译日志存到输出目录
echo "[BUILD] 开始编译 $(date --iso-8601=seconds)"
make -j"${BUILD_JOBS:-6}" > "$OUTDIR/build_stdout.log" 2>&1
BUILD_EXIT=$?
if [ "$BUILD_EXIT" -ne 0 ]; then
    echo "[FATAL] 编译失败，exit=$BUILD_EXIT；详见 $OUTDIR/build_stdout.log" >&2
    echo "$BUILD_EXIT" > "$OUTDIR/exit_code.txt"; exit "$BUILD_EXIT"
fi
echo "[BUILD] 编译成功 $(date --iso-8601=seconds)"

# ====================== ⑥ 运行 ======================
# 命令行参数顺序与单相 poreCementDrive 完全一致（主程序注释里有逐项说明）：
#   N uMax rho_f mu_f v_inf maxT outDir/ nuScale dt_dem totalPV rampTime openOutletPV
#   calibrationMode packFile kMinMd kMaxMd phiMin phiMax
# 两相、润湿、入口控制、驱动协议等新参数不走命令行，C++ 直接从 CONFIG_FILE 读取。
CMD=(mpirun -np "$NP" ./oilwaterDrive
    "$N" "$UMAX" "$RHO_F" "$MU_F" "$V_INF" "$MAX_T" "$OUTDIR/"
    "$NU_SCALE" "$DT_DEM" "$TOTAL_PV" "$RAMP_TIME" "$OPEN_OUTLET_PV" "$CALIBRATION_MODE" "$PACK_FILE"
    "$K_MIN_MD" "$K_MAX_MD" "$PHI_MIN" "$PHI_MAX")
# 把完整命令存成 run_command.sh，以后可以手工复现（%q 会正确转义空格等字符）
printf '%q ' "${CMD[@]}" > "$OUTDIR/run_command.sh"; printf '\n' >> "$OUTDIR/run_command.sh"

# ★★ 2026-10-27：程序的输出全部写进 run_stdout.log，屏幕上原来什么都不显示，看起来像"卡住了"。
#    现在开一个后台小循环，每 60 秒把日志里最新的一条 [step=…] 进度行打印到屏幕上；
#    主程序仍在前台运行，按 Ctrl+C 的效果和以前一样（会中止计算）。
echo "[RUN] 开始计算 $(date --iso-8601=seconds)，每分钟在屏幕上显示一次进度（[step=步数] … PV=已注入量）。"
echo "      完整日志：$OUTDIR/run_stdout.log    ★ 计算过程中不要按 Ctrl+C，也不要关闭这个终端。"
(
    while sleep 60; do
        kill -0 $$ 2>/dev/null || exit 0                                    # 主脚本已结束 → 退出
        LAST=$(grep -E '^\[step=' "$OUTDIR/run_stdout.log" 2>/dev/null | tail -1)
        [ -n "$LAST" ] && echo "[进度 $(date +%H:%M)] $LAST"
        grep -qE '^\[(ERROR|FATAL)\]' "$OUTDIR/run_stdout.log" 2>/dev/null && { echo "[进度] 日志里出现了 [ERROR]/[FATAL]，程序即将停止"; exit 0; }
    done
) &
PROGRESS_PID=$!
"${CMD[@]}" > "$OUTDIR/run_stdout.log" 2>&1           # 真正运行；程序输出全部进 run_stdout.log
PROGRAM_EXIT_CODE=$?
kill "$PROGRESS_PID" 2>/dev/null; wait "$PROGRESS_PID" 2>/dev/null
echo "$PROGRAM_EXIT_CODE" > "$OUTDIR/program_exit_code.txt"
echo "$PROGRAM_EXIT_CODE" > "$OUTDIR/exit_code.txt"

# ====================== ⑦ 汇总 ======================
# 从 run_stdout.log 抓出错误行、完成标志、渗透率结果和最后 20 条 [PV] 进度行，
# 写成 validation_summary.txt，不必翻整个日志就能看出跑没跑完
# ★★ 2026-10-25：旧版只抓 [ERROR]，而启动阶段的退出都是 [FATAL]/[PARAM_ERROR]，
#    所以"错误 7"在 summary 里什么都看不到。现在一起抓，并在屏幕（slurm 输出）上直接打印。
{
    echo "=== errors_and_completion ==="
    grep -E '^\[(ERROR|FATAL|FATAL_GEOMETRY|PARAM_ERROR)\]|^        |达到总注入|本次运行不完整|PERMEABILITY_RESULT' "$OUTDIR/run_stdout.log" | tail -40 || true
    echo "=== preflight ==="
    grep -E '^\[(PRESSURE_PREFLIGHT|TIME_BUDGET|TAU)\]' "$OUTDIR/run_stdout.log" | tail -10 || true
    echo "=== pv_tail ==="
    grep '\[PV\]' "$OUTDIR/run_stdout.log" | tail -20 || true
} > "$OUTDIR/validation_summary.txt"

echo "======================================================"
if [ "$PROGRAM_EXIT_CODE" -eq 0 ]; then
    echo " 运行完成(exit 0)  $(date --iso-8601=seconds)"
else
    echo " 运行异常退出 exit=$PROGRAM_EXIT_CODE  $(date --iso-8601=seconds)"
    echo " ---- 退出原因（摘自 run_stdout.log）----"
    grep -E '^\[(ERROR|FATAL|FATAL_GEOMETRY|PARAM_ERROR)\]|^        ' "$OUTDIR/run_stdout.log" | tail -20 || true
fi
echo " 日志与结果已保存到: $OUTDIR"
echo "======================================================"

# 脚本退出码 = 主程序退出码，方便批处理/作业系统判断成败
exit "$PROGRAM_EXIT_CODE"

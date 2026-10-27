#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_production.py —— 出砂模拟的颗粒级后处理：逐颗追踪弱颗粒的运移、筛前堵塞与产出。

【这个脚本做什么】
  读取 LIGGGHTS 每隔一段时间写出的颗粒快照（post/particles_*.liggghts），
  对每一颗“可动弱颗粒”(type 2) 回放它的整条轨迹，判断它最终是
  已产出 / 运移中 / 滞留 / 筛前堵塞 / 已穿网 / 未显著移动 / 异常丢失，
  并统计产出数随时间的曲线和颗粒数守恒。

【什么时候、由谁调用】
  - 自动：C++ 主程序 oilwaterDrive 在生产模式跑完、dump 转成 CSV 之后，
    会执行 `python3 analyze_production.py <outDir>`；本脚本返回非零
    （异常丢失或颗粒数不守恒）时主程序以退出码 6 结束。
  - 手工：
      单算例：   python3 analyze_production.py <runDir>
      多算例对比：python3 analyze_production.py --compare <runDir1> <runDir2> ...

【输入文件】
  <runDir>/post/<前缀>particles_*.liggghts   颗粒快照（含 particles_final）
  <runDir>/post/<前缀>production_events_p*.csv  主程序在出口删除颗粒前写的“产出事件”
  <runDir>/input_parameters.env / run_stdout.log  取 OUTPUT_PREFIX 和 DEM 步长
  config/defaults.env（经 parameter_config 读取）  取阈值 POST_V_MOVE_M_S 等

【输出文件（单算例，全部写到 <runDir>/analysis/）】
  production_summary.csv     每颗粒一行：首次显著移动/到筛前/堵塞/穿网/到出口/产出时间、
                             最大 x 位移、重启动次数、筛前重启动次数、最终状态
  screen_particle_events.csv 只含到达过筛前统计区的颗粒
  production_curve.csv       t, N_out(t), N_remaining_mobile(t)
  screen_summary.csv         筛前到达/堵塞/重启动/穿网/产出的颗粒数汇总
  production_curves.png      产出数-时间、剩余可动数-时间（无 matplotlib 时跳过）
  summary.txt                文字汇总 + 颗粒数守恒核算
  多算例对比：当前目录 compare_production.png

【状态定义（阈值来自 config：V_MOVE=POST_V_MOVE_M_S，V_STALL=POST_V_STALL_M_S，
  STALL_FRAMES=POST_STALL_FRAMES）】
  已产出     — 主程序在出口删除该颗粒前写出了显式 production event（唯一的“真产出”依据）
  异常丢失   — 末帧里找不到它，却没有产出事件（颗粒从非出口边界消失，必须为 0）
  筛前堵塞   — 曾在筛前统计区内连续 STALL_FRAMES 帧 |v|<V_STALL，且末帧仍停在筛前
  已穿网(未产出) — 球体整个越过筛网平面，但还没走到出口
  未显著移动 — 总位移 < 0.1r；注意：不能据此断言它的胶结键没断
  运移 / 运移(曾重启动) — 末帧仍在动(|v|>V_MOVE)；后者表示中途停过又动起来
  滞留       — 曾显著移动(位移>0.1r)，但末帧 |v|<V_STALL

【关于“重启动”的说明】
  这里的“重启动”(n_restart / n_screen_restart) 只是一个后处理标签：指某颗粒在
  显著移动之后曾经停滞（连续 STALL_FRAMES 帧 |v|<V_STALL），后来又再次运动
  (|v|>V_MOVE)。它描述的是单颗颗粒的运动学历史，与注入方式无关。
  真正的“停注—重启”驱动协议由主程序的 DRIVE_PROTOCOL=stop_restart 实现
  （先注到 STOP_AT_PV，停注保持 HOLD_TIME_S，再恢复注入），其阶段和流量历史
  写在 <runDir>/restart_protocol.csv，不由本脚本处理。

【退出码】
  0 = 颗粒数守恒且异常丢失为 0；4 = 不守恒或有异常丢失；其他异常直接抛出。
"""
import sys, glob, re, csv, os   # sys=命令行参数; glob=按通配符找文件; re=正则; csv=写表; os=路径
from parameter_config import load_parameters, real, text  # 统一从 config/defaults.env 读参数

# ── 后处理阈值：全部来自配置文件，脚本里不写死 ──
CONFIG = load_parameters()
OUTLET_X = real(CONFIG, "DOMAIN_LX_M")                 # 出口截面 x 坐标 = 计算域长度 [m]
SCREEN_X = real(CONFIG, "SCREEN_X_M")                  # 筛网所在 x 坐标 [m]
SCREEN_CAPTURE_BAND = real(CONFIG, "SCREEN_CAPTURE_BAND_M")  # 筛网上游“筛前统计区”的厚度 [m]
_dt_fallback = text(CONFIG, "POST_DT_DEM_FALLBACK_S")  # 日志里找不到步长时的兜底值("auto"=用 DT_DEM)
DT_DEM_DEF = real(CONFIG, "DT_DEM") if _dt_fallback == "auto" else float(_dt_fallback)
V_MOVE = real(CONFIG, "POST_V_MOVE_M_S")               # 速度高于它 → 认为“在运动” [m/s]
V_STALL = real(CONFIG, "POST_V_STALL_M_S")             # 速度低于它 → 认为“停住了” [m/s]
STALL_FRAMES = int(text(CONFIG, "POST_STALL_FRAMES"))  # 连续多少帧低于 V_STALL 才算“停滞”
if (DT_DEM_DEF <= 0 or V_MOVE <= 0 or V_STALL < 0 or V_STALL >= V_MOVE
        or STALL_FRAMES <= 0 or not 0 < SCREEN_X < OUTLET_X or SCREEN_CAPTURE_BAND <= 0):
    raise ValueError("invalid post-processing parameters from run_production.sh")  # 参数不合理就立即停，避免错分类

def parse_dt(run_dir):
    """从运行日志里解析真实的 DEM 步长 t_step（秒/步）；找不到则返回配置缺省值。"""
    prefix = experiment_output_prefix(run_dir)
    candidates = ["run_stdout.log"]
    if prefix:
        candidates.append(os.path.join("post", f"{prefix}log.liggghts"))
    candidates.append("post/log.liggghts")  # 兼容旧版没有前缀的运行目录
    for cand in candidates:  # 依次尝试：主程序 stdout → 带前缀的 LIGGGHTS 日志 → 旧日志
        p = os.path.join(run_dir, cand)
        if os.path.exists(p):
            with open(p, errors="replace") as f:
                for line in f:                            # 找 LIGGGHTS 脚本里的 "variable t_step equal ..."
                    m = re.search(r"variable t_step equal ([0-9.eE+-]+)", line)
                    if m: return float(m.group(1))
    return DT_DEM_DEF                                      # 都没找到 → 用缺省步长

def experiment_output_prefix(run_dir):
    """返回本次运行的输出文件名前缀（OUTPUT_PREFIX，形如 "实验名__"），没有则返回空串。"""
    parameters_path = os.path.join(run_dir, "input_parameters.env")
    if os.path.exists(parameters_path):
        with open(parameters_path, errors="replace") as handle:
            for line in handle:
                if line.startswith("OUTPUT_PREFIX="):
                    return line.rstrip("\n").split("=", 1)[1]
    return os.environ.get("OUTPUT_PREFIX", "")


def load_frames(run_dir):
    """读 post/ 下所有 particles_*.liggghts，按 DEM 步数排序，返回 [(步数, {颗粒id: 数据}), ...]。

    每个颗粒的数据是元组 (x, y, z, vx, vy, vz, radius, type)。
    """
    prefix = experiment_output_prefix(run_dir)
    files = glob.glob(os.path.join(run_dir, "post", f"{prefix}particles_*.liggghts"))
    def dump_step(fn):
        """读 dump 头部第 2 行（ITEM: TIMESTEP 的下一行）得到该帧的步数。"""
        with open(fn, errors="replace") as fh:
            head = [fh.readline() for _ in range(2)]
        return int(head[1].strip())
    files = sorted(files, key=dump_step)                 # 按真实步数排序，这样 particles_final 也能排对位置
    frames = []
    for fn in files:                                      # 一个文件 = 一帧
        step = dump_step(fn)
        atoms = {}                                        # 该帧 {颗粒id: 数据元组}
        with open(fn, errors="replace") as fh:
            lines = fh.readlines()
        for ln in lines[9:]:                              # 前 9 行是 dump 头部，第 10 行起是颗粒数据
            p = ln.split()                                # 列顺序: id type x y z vx vy vz radius ...
            if len(p) >= 9:                               # 字段够多才算有效行（防文件末尾残行）
                atoms[int(p[0])] = (float(p[2]), float(p[3]), float(p[4]),    # x y z
                                    float(p[5]), float(p[6]), float(p[7]),    # vx vy vz
                                    float(p[8]), int(p[1]))                   # radius, type
        frames.append((step, atoms))
    return frames

def load_production_events(run_dir):
    """读取主程序在 delete_atoms 之前写出的显式出口事件，返回 {颗粒id: 产出时刻(s)}。

    每个 MPI 进程各写一个 production_events_p*.csv；同一颗粒若出现多次取最早时刻。
    """
    events = {}
    prefix = experiment_output_prefix(run_dir)
    for fn in glob.glob(os.path.join(run_dir, "post", f"{prefix}production_events_p*.csv")):
        with open(fn, newline="") as fh:
            for row in csv.DictReader(fh):
                pid = int(row["id"])
                events[pid] = min(events.get(pid, float("inf")), float(row["time_s"]))
    return events

def analyze(run_dir):
    """对单个运行目录做完整的颗粒级产出统计，写出 CSV/图/汇总；守恒则返回 0，否则返回 4。"""
    dt = parse_dt(run_dir)                                 # DEM 步长：步数 → 物理时间
    frames = load_frames(run_dir)                          # 所有颗粒快照
    production_events = load_production_events(run_dir)    # 显式产出事件
    if not frames:
        sys.exit(f"未找到 {run_dir}/post/particles_*.liggghts")
    last_step = frames[-1][0]                              # 最后一帧步数，用来判断颗粒是否“中途消失”
    first = frames[0][1]                                   # 第一帧 = 初始状态

    # 可动颗粒 = 弱颗粒(type 2)。骨架(type 1)和筛网(type 3)是冻结的，不参与统计。
    # dump 里没有逐键的断裂信息，所以这里只统计“看得见”的位移、速度和显式产出事件。
    mobile = [pid for pid, d in sorted(first.items()) if d[7] == 2]

    rows, produced_t = [], []                             # rows=每颗一行的统计; produced_t=所有产出时刻
    for pid in mobile:                                    # 逐颗回放它的轨迹
        hist = [(s, a[pid]) for s, a in frames if pid in a]   # 该颗出现过的所有帧 [(步数, 数据)]
        x0, r = hist[0][1][0], hist[0][1][6]              # 初始 x 位置、半径
        t_mobility = t_front = None                        # 首次显著移动时刻 / 首次到达出口截面时刻
        t_screen = t_screen_blocked = t_screen_pass = None  # 首次进入筛前区 / 首次判定筛前堵塞 / 首次穿网
        maxdx, n_restart, n_screen_restart = 0.0, 0, 0     # 最大 x 位移 / 重启动次数 / 筛前重启动次数
        stalled, stall_run = False, 0                      # 停滞状态机：是否已停滞 / 当前连续慢速帧数
        screen_stalled, screen_stall_run = False, 0        # 同上，但只在筛前统计区内计数

        # ── 逐帧扫描：下面是“显著移动 → 停滞 → 重启动”状态机 ──
        # 物理含义：颗粒被水流拖离骨架后可能卡在孔喉里停住(停滞)，之后又被冲开继续走(重启动)。
        for s, d in hist:
            t = s*dt                                       # 步数 → 时间 [s]
            dx = d[0]-x0                                   # 相对初始位置的 x 位移（主流方向）
            maxdx = max(maxdx, dx)
            speed = (d[3]**2+d[4]**2+d[5]**2)**0.5         # 速度大小 |v|
            # (1) 位移超过 0.1 倍半径才算“显著移动”；这只定义运移起点，不代表胶结键断裂时刻
            if t_mobility is None and abs(dx) > 0.1*r: t_mobility = t
            # (2) 球前缘碰到出口截面 → 记录到达出口时刻
            if t_front is None and d[0]+r >= OUTLET_X: t_front = t
            # (3) 筛前统计区：球心在筛网上游、且球前缘已进入筛网上游 SCREEN_CAPTURE_BAND 厚度的带内
            near_screen = d[0] < SCREEN_X and d[0]+r >= SCREEN_X-SCREEN_CAPTURE_BAND
            if t_screen is None and near_screen: t_screen = t
            # (4) 球后缘越过筛网平面 → 整颗已穿网
            if t_screen_pass is None and d[0]-r > SCREEN_X: t_screen_pass = t
            if t_mobility is not None:                     # 只在显著移动之后才统计停滞/重启动
                # 全域停滞计数
                if speed < V_STALL:                        # 这一帧很慢
                    stall_run += 1                         # 连续慢速帧 +1
                    if stall_run == STALL_FRAMES: stalled = True   # 连续够 STALL_FRAMES 帧 → 判定“停滞”
                else:                                      # 这一帧在动
                    if stalled and speed > V_MOVE:         # 之前判定停滞过，现在又超过 V_MOVE → 一次重启动
                        n_restart += 1; stalled = False
                    stall_run = 0                          # 只要动了就重置连续计数
                # 筛前停滞计数：逻辑相同，但要求颗粒位于筛前统计区内（用于识别“筛网堵塞”）
                if near_screen and speed < V_STALL:
                    screen_stall_run += 1
                    if screen_stall_run == STALL_FRAMES:
                        screen_stalled = True
                        if t_screen_blocked is None: t_screen_blocked = t   # 记录首次筛前堵塞时刻
                else:
                    if screen_stalled and speed > V_MOVE:  # 堵塞后又被冲开 → 筛前重启动
                        n_screen_restart += 1
                        screen_stalled = False
                    screen_stall_run = 0

        # ── 用末帧状态给颗粒贴最终标签（判断顺序即优先级）──
        s_end, d_end = hist[-1]                            # 该颗最后一次出现的帧
        v_end = (d_end[3]**2+d_end[4]**2+d_end[5]**2)**0.5 # 末速度
        if pid in production_events:                       # 有显式出口事件 → 真实产出
            state, t_prod = "已产出", production_events[pid]
            produced_t.append(t_prod)
        elif s_end != last_step:                           # 末帧没它、又没有产出事件 → 从非出口边界消失
            state, t_prod = "异常丢失", ""
        elif t_screen_blocked is not None and d_end[0] < SCREEN_X and d_end[0]+r >= SCREEN_X-SCREEN_CAPTURE_BAND and v_end < V_STALL:
            state, t_prod = "筛前堵塞", ""                  # 堵塞过且末帧仍停在筛前
        elif t_screen_pass is not None:
            state, t_prod = "已穿网(未产出)", ""
        elif t_mobility is None:
            state, t_prod = "未显著移动", ""
        elif v_end > V_MOVE:                               # 末帧还在动
            state = "运移(曾重启动)" if n_restart > 0 else "运移"
            t_prod = ""
        else:                                              # 显著移动过但末帧停住
            state, t_prod = "滞留", ""
        rows.append([pid,                                  # 汇总该颗的一行（列顺序与下面表头一致）
                     f"{t_mobility:.4f}" if t_mobility is not None else "",
                     f"{t_screen:.4f}" if t_screen is not None else "",
                     f"{t_screen_blocked:.4f}" if t_screen_blocked is not None else "",
                     f"{t_screen_pass:.4f}" if t_screen_pass is not None else "",
                     f"{t_front:.4f}" if t_front is not None else "",
                     f"{t_prod:.4f}" if t_prod != "" else "",
                     f"{maxdx:.6f}", n_restart, n_screen_restart, state])

    adir = os.path.join(run_dir, "analysis")              # 输出子目录 analysis/
    os.makedirs(adir, exist_ok=True)

    # ── 输出 1：每颗粒明细 ──
    with open(os.path.join(adir, "production_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id","t_mobility_onset_s","t_reach_screen_zone_s","t_screen_blocked_s",
                    "t_pass_screen_s","t_reach_outlet_s","t_produced_s","max_x_disp_m",
                    "n_restart","n_screen_restart","final_state"])
        w.writerows(rows)

    # ── 输出 2：只保留到达过筛前统计区的颗粒 ──
    screen_rows = [row for row in rows if row[2] != ""]
    with open(os.path.join(adir, "screen_particle_events.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id","t_reach_screen_zone_s","t_screen_blocked_s","t_pass_screen_s",
                    "n_screen_restart","t_produced_s","final_state"])
        for row in screen_rows:
            w.writerow([row[0], row[2], row[3], row[4], row[9], row[6], row[10]])

    # ── 输出 3：产出曲线 N_out(t) 与剩余可动颗粒数（阶梯曲线，每产出一颗 +1）──
    produced_t.sort()
    t_end = last_step*dt                                  # 模拟结束时间
    with open(os.path.join(adir, "production_curve.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s","N_out","N_remaining_mobile"])
        w.writerow([0.0, 0, len(mobile)])                 # 起点: t=0 时 0 产出、全部可动
        for i, t in enumerate(produced_t):
            w.writerow([f"{t:.4f}", i+1, len(mobile)-(i+1)])
        w.writerow([f"{t_end:.4f}", len(produced_t), len(mobile)-len(produced_t)])  # 终点

    # ── 输出 4：筛网相关计数汇总 ──
    n_out, n_mob = len(produced_t), len(mobile)           # 产出数 / 可动总数
    states = [r[-1] for r in rows]                        # 所有颗粒的最终状态
    n_screen_arrived = len(screen_rows)                   # 到过筛前统计区
    n_screen_blocked_ever = sum(row[3] != "" for row in rows)      # 曾判定筛前堵塞
    n_screen_passed = sum(row[4] != "" for row in rows)            # 曾穿网
    n_screen_restarted = sum(int(row[9]) > 0 for row in rows)      # 堵塞后曾重启动
    with open(os.path.join(adir, "screen_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mobile_total","reached_screen_zone","blocked_ever","restarted_after_screen_blockage",
                    "passed_screen","produced","final_screen_blocked"])
        w.writerow([n_mob, n_screen_arrived, n_screen_blocked_ever, n_screen_restarted,
                    n_screen_passed, n_out, states.count("筛前堵塞")])

    # ── 输出 5：文字汇总 ──
    summary = (f"运行目录          : {run_dir}\n"
               f"可动颗粒总数      : {n_mob}\n"
               f"已产出            : {n_out}  (总产出率 {100.0*n_out/max(n_mob,1):.1f}%)\n"
               f"运移中            : {states.count('运移')+states.count('运移(曾重启动)')}\n"
               f"  其中曾重启动    : {states.count('运移(曾重启动)')}\n"
               f"到达筛前统计区    : {n_screen_arrived}\n"
               f"曾形成筛前堵塞    : {n_screen_blocked_ever}\n"
               f"堵塞后重启动      : {n_screen_restarted}\n"
               f"已经穿过筛网      : {n_screen_passed}\n"
               f"末态筛前堵塞      : {states.count('筛前堵塞')}\n"
               f"滞留              : {states.count('滞留')}\n"
               f"未显著移动        : {states.count('未显著移动')}\n"
               f"异常丢失          : {states.count('异常丢失')}\n")
    if produced_t:
        summary += f"首颗产出时间      : {produced_t[0]:.3f} s\n"
        summary += f"末颗产出时间      : {produced_t[-1]:.3f} s\n"

    # ── 颗粒数守恒核算：初始可动 = 末帧仍在域内 + 显式产出 + 异常丢失 ──
    # y/z 四个侧面和入口面都是硬反弹边界，颗粒只能从右端出口离开；
    # 所以“异常丢失”理论上必须为 0，否则说明颗粒穿墙或被错误删除。
    OUTLET_HALF = real(CONFIG, "POST_OUTLET_HALF_M")          # 上/下游分界 x（保留参数，供扩展分类用）
    n_remaining = sum(1 for pid in mobile if pid in frames[-1][1])  # 末帧仍在域内的可动颗粒数
    n_abnormal = states.count("异常丢失")
    summary += ("\n── 颗粒数守恒核算 ──────────────\n"
                f"初始可动颗粒(type2) : {n_mob}\n"
                f"域内剩余(末帧)      : {n_remaining}\n"
                f"显式出口产出事件    : {n_out}\n"
                f"异常丢失            : {n_abnormal}\n"
                f"守恒                 : {n_remaining}+{n_out}+{n_abnormal}={n_remaining+n_out+n_abnormal}"
                f" / {n_mob} {'✓' if n_remaining+n_out+n_abnormal == n_mob and n_abnormal == 0 else '⚠'}\n")

    with open(os.path.join(adir, "summary.txt"), "w") as f:
        f.write(summary)
    print(summary)                                        # 同时打印到屏幕（会进 run_stdout.log）

    # ── 输出 6：曲线图（matplotlib 缺失或字体问题时只警告，不影响 CSV）──
    try:
        import matplotlib
        matplotlib.use("Agg")                             # 无界面后端，服务器上也能出图
        import matplotlib.pyplot as plt
        plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC","WenQuanYi Micro Hei","SimHei","DejaVu Sans"]  # 中文字体候选
        plt.rcParams["axes.unicode_minus"] = False        # 正常显示负号
        ts = [0.0]+produced_t+[t_end]                     # 横轴: 0 → 各产出时刻 → 结束
        ns = [0]+list(range(1, n_out+1))+[n_out]          # 纵轴: 累计产出数(阶梯)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))     # 一行两个子图
        ax[0].step(ts, ns, where="post"); ax[0].set_xlabel("t [s]")   # 左: 产出数随时间
        ax[0].set_ylabel("N_out"); ax[0].set_title("Produced particles vs time")
        ax[0].grid(alpha=.3)
        ax[1].step(ts, [n_mob-n for n in ns], where="post", color="tab:red")  # 右: 剩余可动数
        ax[1].set_xlabel("t [s]"); ax[1].set_ylabel("N_remaining")
        ax[1].set_title("Remaining mobile particles vs time"); ax[1].grid(alpha=.3)
        fig.tight_layout()
        fig.savefig(os.path.join(adir, "production_curves.png"), dpi=150)
        print(f"图已保存: {adir}/production_curves.png")
    except Exception as e:
        print(f"[WARN] 绘图失败({e}), CSV 数据仍可用")
    conserved = n_remaining + n_out + n_abnormal == n_mob
    return 0 if conserved and n_abnormal == 0 else 4      # 非零会让主程序以退出码 6 结束

def compare(dirs):
    """多算例对比：把若干运行目录的产出曲线和最终产出率画到同一张图（compare_production.png）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC","WenQuanYi Micro Hei","SimHei","DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    labels, rates, curves = [], [], []                    # 工况名 / 最终产出率 / 曲线数据
    for d in dirs:
        c = os.path.join(d, "analysis", "production_curve.csv")
        if not os.path.exists(c):                          # 还没做过单算例分析 → 先分析一遍
            result = analyze(d)
            if result != 0:
                raise RuntimeError(f"run validation failed for {d}")
        ts, ns, nmob = [], [], 0
        with open(c) as f:
            for row in csv.DictReader(f):
                ts.append(float(row["t_s"])); ns.append(int(row["N_out"]))
                nmob = int(row["N_out"]) + int(row["N_remaining_mobile"])  # 总可动 = 产出 + 剩余
        labels.append(os.path.basename(os.path.normpath(d)))  # 目录名作图例
        rates.append(100.0*ns[-1]/max(nmob,1))            # 最终产出率 %
        curves.append((ts, ns))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for (ts, ns), lb in zip(curves, labels):              # 左图: 各工况产出曲线叠加
        ax[0].step(ts, ns, where="post", label=lb)
    ax[0].set_xlabel("t [s]"); ax[0].set_ylabel("N_out")
    ax[0].set_title("Production curves by case"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].bar(labels, rates, color="tab:blue")            # 右图: 各工况产出率柱状
    ax[1].set_ylabel("Production ratio [%]"); ax[1].set_title("Production ratio by case")
    for i, v in enumerate(rates):                          # 柱顶标数值
        ax[1].text(i, v+1, f"{v:.0f}%", ha="center")
    fig.tight_layout()
    fig.savefig("compare_production.png", dpi=150)        # 保存到当前目录
    print("对比图已保存: compare_production.png")

if __name__ == "__main__":                                # 命令行入口
    if len(sys.argv) < 2:                                  # 没给参数 → 打印用法
        sys.exit(__doc__)
    if sys.argv[1] == "--compare":                        # 多算例对比
        compare(sys.argv[2:])
    else:                                                  # 单个运行目录
        raise SystemExit(analyze(sys.argv[1].rstrip("/")))  # 异常丢失或不守恒时返回非零

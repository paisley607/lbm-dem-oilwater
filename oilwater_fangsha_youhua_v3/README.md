# oilwater_fangsha_youhua — 油水两相 LBM-DEM 出砂模型（审计修订版，2026-10-27 第四次修订）

> **2026-10-27 更新（最新，先读这里）**：本机测试剖面 `config/local_test.env` 已在云端**完整跑通**（第 2433 步注满
> 0.0002 PV，exit 0；孔隙内总密度全程 0.57–1.02，细粒受力无超限，1022 根胶结键中 85 根被毛细力拉断，396 个细粒全部守恒）。
> 最后一个发散点（第 1343 步，骨架颗粒表面的痕量水密度变负）由内核"修订八"解决，日志新增 `[REPAIR]` 行。
> 详见 **`修改说明_20261025.md` 末尾"补充三"**，其中第 15 条是**你在本机 WSL 上一步一步怎么跑**。

> **2026-10-26 更新**：这一版在云端把 Palabos + LIGGGHTS-INL **真正编译并运行**了你的本机测试
> `config/local_test.env`，复现并修掉了"第 7 步受力审计失败"背后的全部原因，详见 **`修改说明_20261025.md` 末尾"补充二"**。要点：
> 1. 侧壁被当成"油水各 1 的流体"、侧壁初始化顺序错（Δρ=−1.37 的来源）、PSM 颗粒算法符号错（上游自带，质量不守恒）——都已修正。
> 2. 润湿性改为文献中的"几何接触角"，**参数直接给角度** `WALL_CONTACT_ANGLE_DEG_TYPE1/2/3`（<90 水湿，>90 油湿）；
>    旧的 `WALL_WATER_FRACTION_*`、`RHO_WALL` 已删除。水湿、油湿现在都稳定。
> 3. 新参数 `DIFFUSION_TAU`（界面动量守恒）、`DIAGNOSTIC_SCAN_STEPS`（排错扫描）；日志新增 `[TIMING]` 各阶段耗时。
> 4. `in.lbdem` 没有改（云端用的是公开版 LIGGGHTS-INL，不认识你那版的 `normalizedCreateDistance`/`pfcDeformability`，
>    测试时临时换成了公开版写法；你那边请继续用原来的 `in.lbdem`）。

> **2026-10-25 更新（先读这里）**：针对超算上的 7 次报错做了第二次修订，详细原因、改动和你的 4 个问题的回答见
> **`修改说明_20261025.md`**。要点：
> 1. 预检/参数检查/MAX_T 检查都改成"说清楚哪个参数、为什么、该改成多少"；`run_production.sh` 会把退出原因打印到 slurm 日志。
> 2. 修了我之前的三个错误：预检用错粘度（应按油）、痕量组分限速让油水界面变成"墙"、M≈20 在旧算法下发散。
>    新增 `MIXTURE_TAU=1`（局部混合松弛时间），`TAU_W_TARGET` 改为 0.7。
> 3. 新增 `plan_run.py`：提交前 1 秒算出能不能过预检、V_INF 最大多少、要算多少步、大概多久。
> 4. **计算量是现在最大的问题**：M≈20、N=40、5 PV 约 5000 万步。请先跑冒烟测试（`submit.sh` 默认就是），
>    再按 `修改说明_20261025.md` 第 5 节选参数。

这是 `oilwater_fangsha` 按《油水两相 LBM-DEM 模型方法学审计》修改后的版本。
每个源文件都加了面向初学者的中文注释；本文件告诉你 **改了什么、为什么改、怎么跑、还有什么没解决**。

> 一句话：原版的两相内核按原参数第 1 步就失稳、毛细力没传给颗粒、润湿性不可控、胶结键强比流体载荷高 4–6 个数量级。
> 本版把这些都改了，并用 Python 小程序验证了新算法；2026-10-26 起 C++ 已在云端（Palabos 2019 公开版 + LIGGGHTS-INL 公开版）
> 完整编译并跑通本机测试剖面。第一次上集群仍请按第 4 节逐步来。

---

## 1. 目录结构

```
oilwater_fangsha_youhua/
├── oilwaterDrive.cpp          # 主程序（两套格子、主循环、入口/出口、协议、诊断）★ 重写
├── oilwaterShanChen3D.h/.hh   # 两相内核：公共速度 + 精确力 + 几何接触角润湿 + 毛细反作用力 ★ 重写
├── diagScan3D.h               # ★★ 2026-10-26 新增：排错用异常格点扫描（DIAGNOSTIC_SCAN_STEPS）
├── oilwaterCoupling3D.h/.hh   # 两相力（PSM 水 + PSM 油 + 毛细反作用）写入 DEM          ★ 重写
├── parameterConfig.h          # C++ 读 config/*.env
├── parameter_config.py        # Python 读 config/*.env（原交付缺失，已补）             ★ 新增
├── in.lbdem / in2.lbdem       # LIGGGHTS 输入脚本（颗粒、接触、平行键、耦合 fix）
├── Makefile.txt（改名为 Makefile 后使用）/ run_production.sh / env.sh
├── config/
│   ├── defaults.env           # 全部参数（★ 两相/协议参数见"油–水两相"一节）
│   ├── production.env / calibration.env / test_1pv.env
│   └── experiments/           # 稀油对照、粗网格探索、润湿对比、重启动（恒流/恒压）、强键对照等剖面
├── data/oilwater_fangsha.porepack   # 颗粒包（未改）
├── lbdem_framework/           # LBDEMcoupling 框架（PSM 耦合，加注释；2026-10-26 修正 PSM 静止方向符号错误）
├── validation/                # ★ numpy 复现新算法的 6 个验证算例（可在任何电脑上跑）
├── plan_run.py                # ★★ 2026-10-25 新增：提交前的可行性计算器
├── submit.sh                  # ★★ 2026-10-25 修订：64 核，默认先跑冒烟测试
├── 修改说明_20261025.md        # ★★ 2026-10-25：报错原因、改动、四个问题的回答
└── *.py                       # 几何校验、键审计、结果分析（逻辑未改，注释重写）
```

## 2. 相对原版的改动（对应审计报告 Phase 4）

| 编号 | 问题（审计结论） | 本版做法 | 文件 |
|---|---|---|---|
| P0-1 | 两相内核无公共速度、注入 τF 而非 F，按原参数首步失稳；τ_w=0.518 不可用 | 公共速度 u′；碰撞前把动量调到 Shan–Doolen 碰撞后的值 j*=(1−ω)j+ωρu′+F；`TAU_W_TARGET`（2026-10-25 起 0.7）自动反算 uMax；★★ 局部混合松弛时间 `MIXTURE_TAU=1`；痕量组分只限制偏离公共速度的部分 | `oilwaterShanChen3D.hh`，`oilwaterDrive.cpp` |
| P0-2 | 格子密度 0.8 产生约 74 kPa 伪压力；粘度比仅 2 | 格子密度油水都取 1；`MU_O=20e-3`（M≈20，可扫描）；启动打印 Bo、Ca、Re、M | `defaults.env`，`oilwaterDrive.cpp` |
| P0-3 | 出口把两相同时设为本体密度（总密度 1.8）；初始界面 1 格阶跃 | 出口只固定总密度，组成取上游一格；`INITIAL_FRONT_WIDTH_CELLS=4` 的 tanh 过渡 | `OilWaterOutletComposition3D`，初始化段 |
| P0-4 | 毛细力没有传给颗粒 | 固体邻格虚拟密度产生的排斥力取反作用累加到颗粒（含力矩），`CAPILLARY_FORCE_ON_PARTICLES=1` | `oilwaterShanChen3D.hh`，`oilwaterCoupling3D.hh` |
| P0-5 | 接触角不可设定，岩心隐式油湿 | `WALL_CONTACT_ANGLE_DEG_TYPE1/2/3`（2026-10-26 起为几何接触角，直接给角度），`validation/v3_contact_angle.py` 检查 | 同上 |
| P0-6 | 键强 60 MPa 比流体载荷高 4–6 个量级；无损伤变量 | 键强改为细粒尺度量级（0.002 MPa，**待实验标定**）；启动打印键破坏力/拖曳/毛细力对比；每个细粒的局部 S_w 与浸水时间写入 `particle_saturation_p*.csv`。**损伤变量本身需改 LIGGGHTS 内核**（见第 6 节） | `defaults.env`，`oilwaterDrive.cpp` |
| P0-7 | 只有恒流量；"重启动"只是后处理标签；细粒限速 0.05 m/s | `INLET_CONTROL=flow|pressure`；`DRIVE_PROTOCOL=stop_restart`（停注—分级重启，写 `restart_protocol.csv`）；限速按格子马赫数自动算 | `OilWaterInlet3D`，`DriveProtocol` |
| P0-8 | 压降/渗透率只用水格子密度与水粘度 | 两相总压 p=c_s²(ρ_w+ρ_o)+c_s²Gρ_wρ_o；分相流量；总流度 + 以 μ_w 折算的表观渗透率 | `OilWaterFields3D`，`perm_porosity.csv` |
| P1 | Smagorinsky 无依据；`parameter_config.py` 缺失；原代码 35 处编译错误 | `SMAGORINSKY_CONSTANT=0`；补 `parameter_config.py`；修复 `Descriptor<T>` 作用域、头文件重复包含、`std::min`、`weakNow` 作用域 | — |

未改动：颗粒包、DEM 接触模型、筛网、PSM 耦合框架本身（审计判定"保留，需论证"）。

## 3. 验证算例（validation/，纯 numpy，笔记本即可跑）

`validation/lbm2p.py` 用 numpy 逐行复现了 C++ 的两相内核（2026-10-25 起默认就是新算法：局部混合松弛时间 +
只限制偏离公共速度；命令前加 `LBM2P_FIXED_TAU=1` 可切回旧算法对比）。先跑它们再上集群：

| 脚本 | 验证什么 | 2026-10-25 结果（G=2.5，新算法） |
|---|---|---|
| `v1_stability_front.py` | 一维油水前缘稳定性、共存密度 | τ 0.6/2.56、0.8/7.925 都稳定；水里溶油 0.15–0.5%，油里溶水 2–4% |
| `v2_laplace.py` | Laplace 定律 → σ_lb，静止液滴稳定性 | τ_w=0.7、0.8（M≈25）R=8–18 格的液滴都稳定；τ_w=0.6 时 R≤11 的液滴发散。液滴会部分溶进油里（见上一行），Δp·R 随半径变化，σ_lb 只能粗估（0.09–0.15），**需要重新标定**（仅影响打印的 Ca） |
| `v3_contact_angle.py` | 输入接触角 → 实测接触角 | 默认 τ 0.7/5.45、WETTING_GAIN=2.4：30→25°，60→58°，90→88°，120→118°（2026-10-26 几何接触角）；★★ 10-27 加入修订八后 45→40.5°，90→88.4°，135→129.0° |
| `v4_layered_poiseuille.py` | 两相分层流 vs 解析解（粘度比、粘性耦合、流量） | 2026-10-26（DIFFUSION_TAU=1）τ 0.7/5.45（粘度比 25）：L2 误差 1.7%，总流量 +0.9%，水层 +0.4%，油层 +3.0%（10-25 版油层 +20.9%）；★★ 10-27 加入修订八后 L2 2.1%，总流量 +0.1%，水层 −0.6%，油层 +3.2%；旧算法 M=1 时界面变成"墙"（误差 96%） |
| `v5_high_viscosity_ratio.py` | 高粘度比下固定 τ 与混合 τ 对比 | 固定 τ：τ_w=0.7 时连粘度比 2.5 都发散；混合 τ：τ_w=0.7 下粘度比 2.5/25 总流量误差 −2.5%/+1.0%，τ_o=19（M≈50）仍稳定 |
| `v6_porous_stability.py` | ★★ 2026-10-26 新增：二维多孔介质（喉道 1–4 格）中的稳定性 | 几何接触角 θ=30°、45°、90°、135°、150° 3000 步全部稳定（10-27 加入修订八后 45°、135° 复测仍稳定）；旧虚拟密度写法 s_wall=0.3/0.7 在 500 步内发散 |

用法：`cd validation && python3 v1_stability_front.py 0.7 5.45 2.5`（参数：τ_w τ_o G）；v2/v3 用环境变量 `TW=0.7 TO=5.45 python3 v3_contact_angle.py`。
换了 G 或 τ 之后应重新跑 v2 并把 σ_lb 写进 `SIGMA_LB_CALIBRATED`（只用于打印 Ca，不影响计算）。

## 4. 在 HPC 上第一次运行（建议顺序）

0. **先把 `Makefile.txt` 改名为 `Makefile`**（写入你电脑时系统不允许直接创建名为 Makefile 的文件，内容完全相同）。
1. `source env.sh`（按你的集群改 module 名和三个路径）。
2. 先确认依赖版本：Palabos-PFM 需含 `SecuredSmagorinskyRegularizedDynamics` 与 `defineDensity`；LIGGGHTS-INL 需含 `cohesion bond` 与 `couple/lb/onetoone`。
3. 编译：`make`。★★ 2026-10-26 起本版已在云端用 Palabos 公开源码（2019）+ LIGGGHTS-INL 公开版完整编译、运行通过
   （`config/local_test.env`）。你的 Palabos-PFM 若是更老的分支，可能出现个别函数签名差异，按报错改即可。
   **先在自己电脑上跑本机测试**：`CONFIG_FILE=config/local_test.env OUTPUT_ROOT=~/LBDEM_runs NP=8 bash run_production.sh`
   （N=20、注 0.0002 PV、约 2430 步；云端 2 核约 1 小时，8 核估计 15–30 分钟），最后应看到"运行完成(exit 0)"。详细步骤见 `修改说明_20261025.md` 补充三第 15 条。
4. 冒烟测试（★★ 2026-10-25 改为 N=20、只注 0.01 PV，约 2–3 小时）：`sbatch submit.sh`（默认就是它），
   或 `CONFIG_FILE=config/test_1pv.env bash run_production.sh`。跑完看日志里的 `MLU/s`。
5. 提交任何正式算例前：`CONFIG_FILE=<剖面> python3 plan_run.py --np <核数> --mlups-per-core <实测MLU/s ÷ 核数>`。
6. 标定：`RUN_MODE=calibration bash run_production.sh`（单相渗透率应落在 3000–4100 mD）。
7. 生产：`bash run_production.sh`；对比实验用 `CONFIG_FILE=config/experiments/xxx.env bash run_production.sh`。
8. 看 `run_stdout.log` 开头的 `═══ 两相参数与无量纲数 ═══`、`[PRESSURE_PREFLIGHT]`、`[TIME_BUDGET]`、`[BOND_SCALE]`，
   运行中看 `[DENSITY]`（实测入口−出口密度差应 ≤0.2）。

## 5. 输出文件（新增/改名的列已在表头写明）

- `dimensionless_numbers.csv`：Δx、Δt、τ、M、G、σ、Ca_model、Ca_real、Bo、Re —— 写方法节直接用。
- `perm_porosity.csv`：按 PV 采样的总流度、表观渗透率（μ_w 折算）、分相流量、核心压降。
- `restart_protocol.csv`：停注/重启各阶段的驱动因子、流量、压降、产出数、运动颗粒数 → 提取重启动阈值。
- `post/particle_saturation_p*.csv`：每个细粒的局部 S_w、累计浸水时间（损伤律的输入）。
- `ibm_force_audit.csv`：新增 `max_capillary_force_N` 列。
- VTK：`velocity_w`、`velocity_o`、`water_saturation`、`pressure_total`（Pa，相对出口）。

## 6. 仍未解决、投稿前必须处理的事

1. **胶结损伤变量（P0-6）**：LIGGGHTS-INL 的 `cohesion bond` 只有固定阈值。要实现"随饱和度/时间演化的键强"，
   需在 LIGGGHTS 源码 `src/cohesion_model_bond.h`（或你版本中的 `pair_gran_*cohesion*`）的键断裂判据处，
   把 `maxSigmaBond` 乘以一个按颗粒存储的折减因子 (1−D)，D 由 `particle_saturation_p*.csv` 同款的 S_w/浸水时间
   通过 `fix property/atom` 传入。本版已把 S_w 和浸水时间算出来了，内核改动约 50 行。
2. **键强标定**：0.002 MPa 只是把量级拉到"流体能冲断"的范围，需用临界流速实验或 Tang 2024 的产出曲线反演。
3. **毛细数与伪速度**：受压差预检限制，现在的流速很低，模型 Ca 约 1e-5～1e-4，已接近真实近井值；但 Shan-Chen 界面旁的
   伪速度（numpy 液滴算例 ~1e-2 格子单位；★★ 2026-10-27 C++ 本机测试：全场平均 1.8e-3、最大约 0.07 格子单位，
   而注入速度只有 7.3e-6）远大于平均孔隙流速，界面附近的流场以伪速度为主。
   界面接触细粒时毛细力比伪速度拖曳大两个量级以上，所以对细粒受力影响有限，但方法节应说明。方法节必须报告 Ca、M、Re、Bo（程序已打印）。
4. **计算量（★★ 2026-10-25）**：每注 1 PV 步数 ≈ 3φ(L²/k)·ν_lb,w·M/Δρ_max，与 N 无关；M≈20 时约 930 万步/PV。
   要做 M≈20 的完整驱替，需要多节点、缩短岩心或改用允许 τ 更接近 0.5 的多相模型（如颜色梯度 + MRT）。
   另外：粘的那一相在壁面旁有表观滑移（流量偏大 10–35%）；油里会溶 2–4% 的水——都是 Shan-Chen 模型的已知局限。
5. **分辨率**：细粒 10 格/直径，审稿要求通常 ≥20；请做 N=40/60/80 的收敛性对比。
6. **Palabos Secured 钳制**：`SecuredSmagorinskyRegularizedDynamics` 会把 |j|>0.2 的动量静默钳制，本版无法计数；
   若 `flow_stability.csv` 的 max_velocity_lb 长期接近 0.2，说明钳制在起作用。
7. **PSM 耦合基准**：单颗粒阻力/升力、双颗粒 DKT 仍需补（审计 Phase 2 第 5、6 项）。
8. **痕量组分密度修复（★★ 2026-10-27 修订八）**：极少数格点的痕量组分密度会被修回 ρ_res（日志 `[REPAIR]`，本机测试约 95 ppm/步），
   相当于很小的人工质量源，方法节一句话说明并报告 ppm 量级。

## 7. 常见问题

- **编译报 `Descriptor` 未声明**：普通函数里要写 `DESCRIPTOR<T>::c`，模板类里才能写 `Descriptor<T>::c`。
- **第 1 步就 `[ERROR] 流场失稳`**：检查 `TAU_W_TARGET` 是否 ≥0.7（0.6 在 numpy 验证中会发散）、G 是否 ≤3.0、`INITIAL_FRONT_WIDTH_CELLS` 是否 >0、`MIXTURE_TAU` 是否为 1。
- **`[FATAL] 预测的入口−出口格子密度差 Δρ=… 超过上限`**：按提示把 V_INF 降到给出的值并把 MAX_T 设成给出的值；或先运行 `python3 plan_run.py`。不要放宽 MAX_PREDICTED_LB_DENSITY_DROP。
- **`[FATAL] MAX_T=… s 不够`**：把 MAX_T 改成提示的值（或减少 TOTAL_PV）。
- **`[PARAM_ERROR] XXX = …`**：每一行就是一个不合法参数和它的允许范围，照着改。
- **`terminate called after throwing … missing parameter XXX`**：你在用旧的 defaults.env，新程序需要 MIXTURE_TAU 等新参数，请用本文件夹的 defaults.env。
- **`[WARN] 键强比流体/毛细载荷高两个数量级以上`**：你在用旧的 60 MPa 键强，流动不会冲断任何键。
- **想退回原版行为对比**：`CAPILLARY_FORCE_ON_PARTICLES=0`、`INITIAL_FRONT_WIDTH_CELLS=0`、`config/experiments/strong_bond_60mpa_reference.env`。

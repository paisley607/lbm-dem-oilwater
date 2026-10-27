# 参数配置说明（oilwater_fangsha_youhua）

## 1. 文件格式

参数文件是最简单的 `KEY=VALUE` 文本（`.env`），每行一个参数，`#` 开头是注释：

```bash
# 入口速度 [m/s]
V_INF=4.5e-5
INLET_CONTROL=flow
```

- `KEY` 只能用大写字母、数字、下划线，且必须以字母开头。
- 值两端的空白和成对引号会被自动去掉（`"abc"` 等价于 `abc`）。
- 不允许写 shell 语法（`$VAR`、`$(cmd)`、`export`），因为加载器**不会 source 这个文件**，
  而是逐行解析。这样 bash（`config/load_parameters.sh`）、C++（`parameterConfig.h`）和
  Python（`parameter_config.py`）三方都能用同一套规则读到完全相同的值。

## 2. 三层优先级与继承

```
环境变量（最高） > 实验文件（及其 BASE_CONFIG 链） > config/defaults.env（最低）
```

1. **`config/defaults.env` 是全部合法参数的清单，也是唯一的默认值来源。**
   实验文件里出现任何 defaults.env 没有的 KEY，加载器立刻报
   `[FATAL] unknown parameter`，防止拼错参数名却静默生效。
2. **实验文件只写需要改动的参数**，放在 `config/experiments/`。它可以用
   `BASE_CONFIG=相对或绝对路径` 继承另一个文件（相对路径相对于本文件所在目录）。
   加载顺序是 defaults → 基础文件（可多级，先父后子）→ 本文件，后加载的覆盖先加载的。
   循环继承（A 继承 B、B 又继承 A）会被检测并报错。
3. **命令行环境变量永远不会被文件覆盖**，适合临时改一个参数做对比：
   ```bash
   CONFIG_FILE=config/experiments/light_oil_M2.env V_INF=3e-4 bash run_production.sh
   ```
   加载完成后，实际的文件链会导出到 `PARAMETER_CONFIG_CHAIN`（冒号分隔）。

典型用法：

```bash
# 用默认参数
bash run_production.sh
# 用某个实验文件
CONFIG_FILE=config/experiments/light_oil_M2.env bash run_production.sh
# 提交前先估算（★★ 2026-10-25 新增）
CONFIG_FILE=config/experiments/light_oil_M2.env python3 plan_run.py
# 标定初始渗透率
RUN_MODE=calibration bash run_production.sh
```

## 3. 单位后缀约定

单位直接写在参数名后缀里，避免看错量级：

| 后缀 | 单位 |
|---|---|
| `_M`, `_M_S`, `_M_S2` | 米、米/秒、米/秒² |
| `_N`, `_NM`, `_PA`, `_PA_S`, `_RAD_S` | 牛、牛·米、帕、帕·秒、弧度/秒 |
| `_UM`, `_UM_S` | 微米、微米/秒 |
| `_MM`, `_MM_S` | 毫米、毫米/秒 |
| `_GPA`, `_MPA`, `_MD` | 吉帕、兆帕、毫达西 |
| `_KG_M3` | 千克/立方米 |
| `_S` | 秒 |
| `_CELLS`, `_LB_STEPS` | 格子数、LBM 步数（无量纲整数） |
| 无后缀（`N`, `NP`, 比值, 孔隙度, PV, 容差） | 无量纲，除非该行注释另有说明 |

## 4. 本版新增的两相 / 驱动协议参数

下面这些是审计修订版新增的参数，详细说明见 `config/defaults.env` 里对应的注释行
（带 ★ 的）。原版单相参数（N、UMAX、V_INF、DEM_*、BOND_* 等）含义不变。

| 参数 | 一句话含义 |
|---|---|
| `TAU_W_TARGET` | 水相格子松弛时间目标值；>0 时主程序反算 uMax 使 τ_w 等于它（★★ 2026-10-25 起默认 0.7；numpy 验证 0.6 的液滴/接触角算例发散），0 = 沿用 UMAX |
| `MIXTURE_TAU` | ★★ 1 = 局部混合松弛时间（按含水比例混合粘度，两相共用 τ；M 较大时必须），0 = 旧做法（水/油各用固定 τ，τ_o≥4.5 发散） |
| `MIX_SAT_LOW` / `MIX_SAT_HIGH` | ★★ 含水比例低于/高于它们时按纯油/纯水粘度，中间调和平均（默认 0.05/0.95） |
| `TRACE_VELOCITY_CAP` | ★★ 痕量组分偏离公共速度的上限（格子单位，默认 0.3），防止只有残余密度的组分速度失控 |
| `MAX_RUNTIME_LB_DENSITY_RISE` | ★★ 运行中实测入口−出口格子密度差的上限（默认 0.35），超过就停算 |
| `MAX_PREDICTED_LB_DENSITY_DROP` | 启动预检：预测的入口−出口格子密度差上限（★★ 默认 0.2，允许最大 0.35；不是失稳报警，见 修改说明_20261025.md 问题 4） |
| `WALL_CONTACT_ANGLE_DEG_TYPE1` | ★★ 2026-10-26 新：type1 骨架颗粒的接触角 [度]，从水相一侧量：<90 水湿，90 中性，>90 油湿；可用 20～160 |
| `WALL_CONTACT_ANGLE_DEG_TYPE2` | 同上，针对 type2 弱胶结细粒 |
| `WALL_CONTACT_ANGLE_DEG_TYPE3` | 同上，针对 type3 出口筛网 |
| `CONTAINER_WALL_CONTACT_ANGLE_DEG` | 计算域四个侧壁（方形外壳）的接触角，默认 90（中性），一般不用改 |
| `WETTING_GAIN` | 几何接触角公式的标定系数，默认 2.4（numpy 液滴算例标定：输入 60° → 实测 58°）；一般不用改 |
| `DIFFUSION_TAU` | 油水两相之间动量交换的松弛时间，默认 1.0；0 = 旧做法（界面动量不守恒，只做对比用） |
| `SHAN_CHEN_SOLID_THRESHOLD` | 固含率超过它的格点在两相作用（Shan-Chen 力、润湿）里当固体，默认 0.5，一般不用改 |
| `PARTICLE_MACH_LIMIT` | 细粒速度上限（格/步），默认 0.05（2026-10-26 由 0.3 改小：PSM 要求颗粒每步移动远小于 1 格） |
| `DIAGNOSTIC_SCAN_STEPS` | 排错开关：前 N 步扫描全场可疑格点并按位置分类打印 [DIAG] 行；正常计算保持 0 |
| `INITIAL_FRONT_WIDTH_CELLS` | 初始油水界面 tanh 平滑过渡的半宽 [格]，默认 2；过渡区总密度按"压力处处相等"自动给出；0 = 一格阶跃 |
| ~~`RHO_WALL`、`WALL_WATER_FRACTION_TYPE1/2/3`、`CONTAINER_WALL_WATER_FRACTION`~~ | 2026-10-26 已删除（旧虚拟密度润湿模型），实验文件里还写着它们会报 unknown parameter，请改用上面的接触角 |
| `CAPILLARY_FORCE_ON_PARTICLES` | 1 = 把毛细反作用力传给颗粒（默认）；0 = 只传动量交换力（等价原版，用于对比） |
| `INLET_CONTROL` | 入口控制方式：`flow` = 恒流量（入口速度 V_INF）；`pressure` = 恒压力（入口压差 INLET_PRESSURE_PA） |
| `INLET_PRESSURE_PA` | 恒压模式下入口相对出口的压差 [Pa]，只在 INLET_CONTROL=pressure 时使用 |
| `DRIVE_PROTOCOL` | 驱动协议：`continuous` = 缓启后持续驱替；`stop_restart` = 停注—分级重启（下面 5 个参数生效） |
| `STOP_AT_PV` | stop_restart 协议：注到多少 PV 时停注 |
| `HOLD_TIME_S` | stop_restart 协议：停注后保持静止的时间 [s] |
| `RESTART_STAGES` | stop_restart 协议：重启时分几级台阶升流量/压差 |
| `RESTART_STAGE_HOLD_S` | stop_restart 协议：每级台阶保持的时间 [s]，每级的流量/压差/产出数写到 restart_protocol.csv |
| `RESTART_MAX_FACTOR` | stop_restart 协议：最后一级升到原流量/压差的多少倍 |
| `PARTICLE_MACH_LIMIT` | 细粒速度上限对应的格子马赫数：上限 = PARTICLE_MACH_LIMIT × Δx/Δt（DEM_MAX_PARTICLE_SPEED_M_S=0 时自动启用） |
| `SIGMA_LB_CALIBRATED` | 由 validation/v2_laplace.py 标定的格子界面张力（旧算法 G=2.5、τ 0.8/2.5 → 约 0.085；★★ 新算法液滴会部分溶入油中，只能粗估 0.09–0.15），只用于启动时打印模型毛细数 Ca_model；换 G 或 τ 后需重新标定 |
| `SIGMA_PHYS_REFERENCE_N_M` | 真实油水界面张力 [N/m]，只用于 Bond 数、Ca_real 和键强量级提示 |

## 5. 加载器会拒绝什么

- 格式不对的行（没有 `=`、KEY 含小写）
- defaults.env 里没登记的 KEY
- 文件不存在、BASE_CONFIG 循环继承

此外 `run_production.sh` 在编译和 mpirun 之前还会做取值范围检查
（`validate_water_fangsha_parameter_values`）和颗粒包几何校验。

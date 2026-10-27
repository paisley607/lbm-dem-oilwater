#!/bin/bash
# ============================================================================
# config/load_parameters.sh —— 安全读取 KEY=VALUE 参数文件（不 source、不执行任何 shell 代码）
#
# 【为什么不直接 source .env】source 会执行文件里的任意命令，而且拼错的参数名会被静默接受。
#   这里改为逐行解析，并实现下面三条规则：
#   1. defaults.env 是"合法参数名清单"：只有在 defaults.env 里出现过的 KEY 才允许出现在
#      实验文件里，拼错立刻报 [FATAL] unknown parameter。
#   2. BASE_CONFIG 继承：实验文件可以写 BASE_CONFIG=../production.env，加载顺序是
#      defaults → 基础文件(递归，先父后子) → 本文件，后加载的覆盖先加载的；检测循环继承。
#   3. 环境变量优先：调用脚本前已经在环境里的变量（例如 NP=8 bash run_production.sh）
#      不会被任何文件覆盖，方便临时改参数而不动文件。
#
# 【对外接口】
#   load_water_fangsha_parameters <defaults.env> [实验.env]
#       读参数并 export 成 shell 变量；同时 export PARAMETER_CONFIG_CHAIN（冒号分隔的加载链）。
#   validate_water_fangsha_parameter_values
#       对关键参数做取值范围检查（正数/非负/正整数/几何顺序）。
# C++ 侧的 parameterConfig.h 用完全相同的规则再读一遍，保证 bash 与主程序看到同一套值。
# ============================================================================

set -u   # 引用未定义变量即报错，避免拼错变量名悄悄变成空串

# 去掉值两端的空白和成对的引号："abc" / 'abc' → abc
trim_parameter_value() {
    local value="$1"
    value="${value#${value%%[![:space:]]*}}"      # 去左侧空白
    value="${value%${value##*[![:space:]]}}"      # 去右侧空白
    if [[ "$value" == \"*\" && "$value" == *\" ]]; then
        value="${value:1:${#value}-2}"
    elif [[ "$value" == \'*\' && "$value" == *\' ]]; then
        value="${value:1:${#value}-2}"
    fi
    printf '%s' "$value"
}

# 列出一个参数文件里出现的全部 KEY（每行一个）；跳过空行和 # 注释
# KEY 必须是大写字母开头、只含大写字母/数字/下划线
parameter_file_keys() {
    local file="$1" line key value
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line#${line%%[![:space:]]*}}"
        [[ -z "$line" || "${line:0:1}" == "#" ]] && continue
        if [[ "$line" =~ ^([A-Z][A-Z0-9_]*)=(.*)$ ]]; then
            printf '%s\n' "${BASH_REMATCH[1]}"
        fi
    done < "$file"
}

# 从参数文件里取某个 KEY 的值（第一次出现的）；找不到返回 1
# 主要用来读 BASE_CONFIG，在正式加载之前先知道要继承谁
parameter_file_value() {
    local file="$1" wanted="$2" line key value
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line#${line%%[![:space:]]*}}"
        [[ -z "$line" || "${line:0:1}" == "#" ]] && continue
        [[ "$line" =~ ^([A-Z][A-Z0-9_]*)=(.*)$ ]] || continue
        key="${BASH_REMATCH[1]}"
        [[ "$key" == "$wanted" ]] || continue
        trim_parameter_value "${BASH_REMATCH[2]}"
        return 0
    done < "$file"
    return 1
}

# 判断某个 KEY 是否在 defaults.env 里登记过（返回 0 = 已知）
parameter_key_is_known() {
    local defaults_file="$1" key="$2" candidate found=1
    [[ "$key" =~ ^[A-Z][A-Z0-9_]*$ ]] || return 1
    while IFS= read -r candidate; do
        if [[ "$candidate" == "$key" ]]; then found=0; break; fi
    done < <(parameter_file_keys "$defaults_file")
    return "$found"
}

# 加载单个参数文件：逐行解析、检查 KEY 合法、赋值并 export
# 依赖两个全局关联数组（由 load_water_fangsha_parameters 准备）：
#   _PARAMETER_ALLOWED[KEY]  = 1 表示 KEY 在 defaults.env 里登记过
#   _PARAMETER_EXTERNAL[KEY] = 1 表示 KEY 在调用前已存在于环境中（不可被文件覆盖）
load_parameter_file() {
    local file="$1" line key value
    [[ -f "$file" ]] || {
        echo "[FATAL] parameter file not found: $file" >&2
        return 1
    }
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line#${line%%[![:space:]]*}}"
        [[ -z "$line" || "${line:0:1}" == "#" ]] && continue
        [[ "$line" =~ ^([A-Z][A-Z0-9_]*)=(.*)$ ]] || {           # 不是 KEY=VALUE 形式就是格式错误
            echo "[FATAL] invalid parameter line in $file: $line" >&2
            return 1
        }
        key="${BASH_REMATCH[1]}"
        if [[ -z "${_PARAMETER_ALLOWED[$key]+x}" ]]; then        # 未登记的 KEY：多半是拼写错误
            echo "[FATAL] unknown parameter '$key' in $file; check spelling against defaults.env" >&2
            return 1
        fi
        value="$(trim_parameter_value "${BASH_REMATCH[2]}")"
        # 环境变量优先：调用前已设置的变量不改；否则文件值（后加载的覆盖先加载的）生效
        if [[ -z "${_PARAMETER_EXTERNAL[$key]+x}" ]]; then
            printf -v "$key" '%s' "$value"                        # 等价于 KEY="$value"，但 KEY 名来自变量
            export "$key"
        fi
    done < "$file"
}

# 主入口：load_water_fangsha_parameters <defaults.env> [实验.env]
# 步骤：登记合法 KEY → 记录已有环境变量 → 加载 defaults → 递归加载实验文件及其 BASE_CONFIG
load_water_fangsha_parameters() {
    local defaults_file="$1" profile_file="${2:-}" line source base canonical
    [[ -f "$defaults_file" ]] || {
        echo "[FATAL] parameter defaults file not found: $defaults_file" >&2
        return 1
    }
    # 合法 KEY 清单 = defaults.env 里的全部 KEY
    declare -gA _PARAMETER_ALLOWED=()
    while IFS= read -r key; do _PARAMETER_ALLOWED[$key]=1; done < <(parameter_file_keys "$defaults_file")
    # 记录哪些 KEY 在调用前已经是环境变量（这些值"冻结"，任何文件都不能覆盖）
    declare -gA _PARAMETER_EXTERNAL=()
    while IFS= read -r key; do
        if declare -p "$key" >/dev/null 2>&1; then
            _PARAMETER_EXTERNAL[$key]=1
        fi
    done < <(parameter_file_keys "$defaults_file")

    # 先加载 defaults（只给环境里没有的变量赋值）
    load_parameter_file "$defaults_file" || return 1
    declare -gA _PARAMETER_SEEN=()            # 已进入加载链的文件（规范路径），用于检测循环继承
    declare -ga _PARAMETER_CONFIG_CHAIN=()    # 实际加载顺序，最后导出为 PARAMETER_CONFIG_CHAIN
    canonical="$(cd "$(dirname "$defaults_file")" && pwd -P)/$(basename "$defaults_file")"
    _PARAMETER_CONFIG_CHAIN+=("$canonical")
    # 递归加载实验文件：先加载它的 BASE_CONFIG（可能再有自己的 BASE_CONFIG），再加载它本身
    load_parameter_profile_recursive() {
        local file="$1" base_path current_canonical
        [[ -f "$file" ]] || {
            echo "[FATAL] parameter profile not found: $file" >&2
            return 1
        }
        current_canonical="$(cd "$(dirname "$file")" && pwd -P)/$(basename "$file")"
        if [[ -n "${_PARAMETER_SEEN[$current_canonical]+x}" ]]; then     # A 继承 B、B 又继承 A → 死循环，直接报错
            echo "[FATAL] cyclic BASE_CONFIG inheritance detected at $current_canonical" >&2
            return 1
        fi
        _PARAMETER_SEEN[$current_canonical]=1
        base_path="$(parameter_file_value "$file" BASE_CONFIG || true)"
        if [[ -n "$base_path" ]]; then
            if [[ "$base_path" != /* ]]; then                   # 相对路径相对于当前实验文件所在目录
                base_path="$(dirname "$file")/$base_path"
            fi
            load_parameter_profile_recursive "$base_path" || return 1
        fi
        load_parameter_file "$file" || return 1                 # 父文件加载完再加载自己，自己的值覆盖父文件
        _PARAMETER_CONFIG_CHAIN+=("$current_canonical")
    }
    # 指定了实验文件且不是 defaults 本身时才做覆盖加载
    if [[ -n "$profile_file" && "$profile_file" != "$defaults_file" ]]; then
        load_parameter_profile_recursive "$profile_file" || return 1
    fi
    # 把加载链拼成 "a.env:b.env:c.env" 导出，写进日志可追溯参数来源
    PARAMETER_CONFIG_CHAIN="$(IFS=:; printf '%s' "${_PARAMETER_CONFIG_CHAIN[*]}")"
    export PARAMETER_CONFIG_CHAIN
    unset _PARAMETER_ALLOWED _PARAMETER_SEEN
    unset _PARAMETER_EXTERNAL
}

# ---------------- 取值检查工具（参数名作为参数传入，用 ${!1} 间接取值） ----------------

# 要求是有限正数（允许科学计数法，如 1.5e-3）
parameter_require_positive() {
    local key="$1" value="${!1-}"
    awk -v value="$value" 'BEGIN {
        if (value !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ || (value+0) <= 0) exit 1
    }' || {
        echo "[FATAL] parameter $key must be a finite positive number, got '$value'" >&2
        return 1
    }
}

# 要求是有限非负数（允许 0）
parameter_require_nonnegative() {
    local key="$1" value="${!1-}"
    awk -v value="$value" 'BEGIN {
        if (value !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ || (value+0) < 0) exit 1
    }' || {
        echo "[FATAL] parameter $key must be a finite non-negative number, got '$value'" >&2
        return 1
    }
}

# 要求是正整数（不允许前导 0、小数、科学计数法）
parameter_require_integer_positive() {
    local key="$1" value="${!1-}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
        echo "[FATAL] parameter $key must be a positive integer, got '$value'" >&2
        return 1
    }
}

# 对关键参数做范围检查；最后一条检查几何顺序：0 ≤ CORE_X0 < CORE_X1 < SCREEN_X < DOMAIN_LX，
# 以及 0 < PHI_MIN < PHI_MAX < 1
validate_water_fangsha_parameter_values() {
    local key
    for key in N NP BUILD_JOBS SCREEN_WIRE_COUNT DEM_NEIGHBOR_EVERY_STEPS DEM_NEIGHBOR_ONE DEM_NEIGHBOR_PAGE DEM_THERMO_EVERY PARTICLE_CLEANUP_INTERVAL_LB_STEPS; do
        parameter_require_integer_positive "$key" || return 1
    done
    for key in RHO_F MU_F V_INF MAX_T NU_SCALE DT_DEM TOTAL_PV RAMP_TIME MIN_WEAK_DIAMETER_UM DOMAIN_LX_M DOMAIN_LY_M DOMAIN_LZ_M CORE_X1_M SCREEN_X_M RHO_O MU_O RHO_RESIDUAL INITIAL_FRONT_POSITION_M; do
        parameter_require_positive "$key" || return 1
    done
    for key in EQUILIBRATION_TIME_S OPEN_OUTLET_PV DEM_NEIGHBOR_SKIN_M DEM_NEIGHBOR_DELAY_STEPS SHAN_CHEN_G; do
        parameter_require_nonnegative "$key" || return 1
    done
    awk -v lo="$PHI_MIN" -v hi="$PHI_MAX" -v x0="$CORE_X0_M" -v x1="$CORE_X1_M" -v sx="$SCREEN_X_M" -v lx="$DOMAIN_LX_M" 'BEGIN {
        if (!(lo > 0 && hi > lo && hi < 1 && x0 >= 0 && x1 > x0 && x1 < sx && sx < lx)) exit 1
    }' || {
        echo "[FATAL] invalid porosity or x-domain ordering: PHI_MIN/PHI_MAX and CORE_X0_M<CORE_X1_M<SCREEN_X_M<DOMAIN_LX_M are required" >&2
        return 1
    }
    return 0
}

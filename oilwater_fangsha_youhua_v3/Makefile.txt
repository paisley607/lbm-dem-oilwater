# ============================================================================
# Makefile —— oilwater_fangsha_youhua 的编译入口
#
# 【它做什么】在目录下敲 `make`（或由 run_production.sh 自动调用）时，把主程序
#   oilwaterDrive.cpp 和自带耦合框架的两个 .cpp 一起交给 Palabos 的 SCons 构建系统，
#   链接 LIGGGHTS 动态库，生成可执行文件 oilwaterDrive。
#
# 【为什么是 SCons 而不是直接 g++】Palabos 官方用 SCons 管理自己几百个源文件的编译选项，
#   项目只需告诉它"我的源文件是谁、头文件在哪、要链哪个库"，其余由 Palabos 的 SConstruct
#   处理。下面的 SConsArgs 就是这些"告诉它"的参数。
#
# 【为什么用自带的 lbdem_framework】原始 LBDEMcoupling 框架在本项目里做过修改
#   （例如给耦合域补充了 setSpheres 重载），为了保证交付目录可独立复现，直接把改过的
#   框架随目录一起带着，不再依赖外部 LBDEM_ROOT。
#
# 【需要的环境变量】（由 run_production.sh 给默认值，或先 source env.sh）
#   PALABOS_ROOT  : Palabos-PFM 根目录（里面有 SConstruct 和 src/）
#   LIGGGHTS_ROOT : LIGGGHTS-INL 根目录（src/ 下有 libliggghts.so）
#   换机器只需在命令行覆盖，例如：
#     PALABOS_ROOT=/opt/palabos LIGGGHTS_ROOT=/opt/liggghts make
# 【可选】BUILD_JOBS : 并行编译任务数（默认 6）
# ============================================================================

# 把配方行的前缀从 Tab 改成 ">"，避免编辑器把 Tab 换成空格导致 make 报错
.RECIPEPREFIX := >

palabosRoot := $(PALABOS_ROOT)
liggghtsRoot := $(LIGGGHTS_ROOT)
BUILD_JOBS ?= 6

# 自带框架目录：头文件 + 两个框架 .cpp（LIGGGHTS 封装、格子分块）
frameworkDir := lbdem_framework

# 参与编译的源文件：主程序 + 两个框架源文件
projectFiles := oilwaterDrive.cpp \
                $(frameworkDir)/liggghtsCouplingWrapper.cpp \
                $(frameworkDir)/latticeDecomposition.cpp

# SCons 命令：-j 并行数，-f 指定 Palabos 的 SConstruct 文件
SCons := scons -j $(BUILD_JOBS) -f $(palabosRoot)/SConstruct

# 传给 Palabos SConstruct 的参数（注意：续行 "\" 之后不能加注释，所以说明写在这里）
#   palabosRoot   Palabos 根目录
#   projectFiles  上面列出的源文件
#   optimize/debug/profile  三选一的构建类型，这里用优化版(-O3)
#   MPIparallel=true        用 mpicxx 编译，生成 MPI 并行程序
#   SMPparallel=false       不用 OpenMP 线程并行
#   usePOSIX=true           允许用 POSIX 接口（realpath 等）
#   compileFlags            附加编译选项；-std=gnu++11 保证 C++11 语法（auto、>> 嵌套模板等）
#                           在老版本 g++ 默认标准下也能编译
#   libraryPaths/libraries  LIGGGHTS 动态库位置和库名（-L ... -lliggghts）
#   includePaths            头文件搜索路径：自带框架目录、交付根目录
#                           （oilwaterShanChen3D.h / oilwaterCoupling3D.h / parameterConfig.h）、LIGGGHTS src
SConsArgs := palabosRoot=$(palabosRoot) \
             projectFiles="$(projectFiles)" \
             optimize=true \
             debug=false \
             profile=false \
             MPIparallel=true \
             SMPparallel=false \
             usePOSIX=true \
             serialCXX=g++ \
             parallelCXX=mpicxx \
             compileFlags="-Wall -Wnon-virtual-dtor -Wno-write-strings -std=gnu++11" \
             linkFlags="" \
             optimFlags="-O3" \
             debugFlags="-g" \
             profileFlags="-pg" \
             libraryPaths="$(liggghtsRoot)/src" \
             includePaths="$(frameworkDir) . $(liggghtsRoot)/src $(palabosRoot)/externalLibraries $(palabosRoot)/externalLibraries/tinyxml" \
             libraries="liggghts"

.PHONY: all compile clean

# `make` 默认目标 = compile
all: compile

# 编译：调用 SCons；SCons 自己判断哪些 .o 需要重编
compile:
>$(SCons) $(SConsArgs)

# 清理：让 SCons 删掉它生成的中间文件，再删掉残留的 .o 和可执行文件
clean:
>$(SCons) -c $(SConsArgs)
>rm -f *.o oilwaterDrive

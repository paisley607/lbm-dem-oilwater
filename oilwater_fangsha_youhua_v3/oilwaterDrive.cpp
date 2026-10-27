/* =====================================================================================
   oilwaterDrive.cpp —— 油水两相 LBM-DEM 出砂模拟主程序（审计修订版 oilwater_fangsha_youhua）
   =====================================================================================

   【这个程序做什么】
     在一个 2.3 mm × 1 mm × 1 mm 的小岩心里（固定骨架颗粒 + 弱胶结细粒 + 出口筛网），
     用格子 Boltzmann 方法（LBM，Palabos 库）算油、水两相流动，用离散元（DEM，LIGGGHTS 库）
     算细粒的脱落、运移、在筛网处桥堵或产出。两者每个时间步交换一次力和位置。

   【程序的骨架】（和原版相同，便于对照）
     1. 读参数（config/defaults.env + 命令行）→ 2. 单位换算 → 3. 建 DEM 体系（in.lbdem）
     → 4. 建两套 LBM 格子（水 latticeW、油 latticeO）→ 5. 边界与初始相分布
     → 6. 主循环：映射颗粒到格子 → 入口/出口条件 → 两相相互作用 → 碰撞迁移
        → 流体力写入 DEM → DEM 子循环 → 诊断输出 → 7. 收尾与后处理。

   【相对审计前版本改了什么】（编号对应审计报告 Phase 4 的 P0/P1 条目）
     P0-1  两相内核：公共速度 + 精确力注入（oilwaterShanChen3D.hh）；τ_w 由参数 TAU_W_TARGET
           决定（默认 0.8，程序自动反算 uMax），不再沿用单相的 τ≈0.518。
     P0-2  格子密度油水都取 1（密度比的影响用 Bond 数论证，启动时打印）；粘度比由 MU_O 设定。
     P0-3  出口只固定总压，各相密度按上游局部组成分配；初始界面用 tanh 平滑过渡。
     P0-4  毛细反作用力传给颗粒（oilwaterCoupling3D.hh）。
     P0-5  润湿性参数（2026-10-26 起为接触角 WALL_CONTACT_ANGLE_DEG_TYPE1/2/3），可按颗粒类型设定。
     P0-6  胶结损伤：LIGGGHTS 的平行键内核没有损伤变量，这里做了两件能做的事——
           ① 每个细粒附近的局部含水饱和度 S_w 和"浸水时间"逐步记录到 CSV，供损伤律使用；
           ② 键强参数改为细粒尺度量级并在启动时与流体/毛细载荷比较（见 defaults.env 说明）。
           真正的"随饱和度演化的键强"需要改 LIGGGHTS 内核，README 里给了接口位置。
     P0-7  入口支持恒流量（flow）和恒压力（pressure）两种控制；支持"停注—分级重启"协议；
           细粒速度上限改为按格子马赫数自动计算。
     P0-8  压降用两相总压 p = c_s²(ρ_w+ρ_o) + c_s²Gρ_wρ_o；流量分相统计；给出总流度与表观渗透率。
     P1    Smagorinsky 常数默认 0；启动时打印 Ca、Re、M、Bo 等无量纲数；原代码的编译错误
           （Descriptor 作用域、头文件重复包含、std::min、weakNow 作用域）已修复。

   【命令行】（与原版完全一致，脚本可复用）
     mpirun -np NP ./oilwaterDrive N uMax rho_f mu_f v_inf maxT outDir/ \
            [nuScale dt_dem totalPV rampTime openOutletPV calibrationMode packFile kMinMd kMaxMd phiMin phiMax]
     两相与协议参数从 config/defaults.env 读取（见该文件"油–水两相"一节）。
   ===================================================================================== */

#define LBDEM_USE_WEIGHING             // 启用固含率加权（PSM/IMB 风格）耦合，见 ibDef.h

#include "palabos3D.h"
#include "palabos3D.hh"
#include "plb_ib.h"
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include "parameterConfig.h"
#include "liggghtsCouplingWrapper.h"
#include "latticeDecomposition.h"
#include "oilwaterShanChen3D.h"        // 两相 Shan-Chen（公共速度 + 润湿 + 毛细反作用力）
#include "oilwaterCoupling3D.h"        // 两相力汇总写入 DEM
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_lb_coupling_onetoone.h"
#include "diagScan3D.h"               // ★★ 2026-10-26：排错用异常格点扫描（DIAGNOSTIC_SCAN_STEPS>0 时启用）

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor     // D3Q19：三维、19 个离散速度方向

static water_fangsha::ParameterConfig* activeParameterConfig = 0;

// ── 参数读取小工具（从 config/defaults.env 或环境变量取值）──────────────────────────
T environmentDouble(char const* name) {
    if (!activeParameterConfig) throw std::runtime_error("parameter configuration is not initialized");
    return (T)activeParameterConfig->real(name);
}
plint environmentInteger(char const* name) {
    if (!activeParameterConfig) throw std::runtime_error("parameter configuration is not initialized");
    return (plint)activeParameterConfig->integer(name);
}
std::vector<T> environmentDoubleList(char const* name) {
    if (!activeParameterConfig) throw std::runtime_error("parameter configuration is not initialized");
    std::vector<double> rawValues = activeParameterConfig->realList(name);
    return std::vector<T>(rawValues.begin(), rawValues.end());
}
std::string environmentText(char const* name) {
    if (!activeParameterConfig) throw std::runtime_error("parameter configuration is not initialized");
    return activeParameterConfig->text(name);
}

// ── 格点有效孔隙率（壁面 0；颗粒内部按 1−固含率）────────────────────────────────────
inline T effectivePorosity(Cell<T,DESCRIPTOR>& cell) {
    if (dynamic_cast<BounceBack<T,DESCRIPTOR>*>(&cell.getDynamics()) != 0) return (T)0.;
    IBdynamicsParticleData<T,DESCRIPTOR>* particle = getParticleDataFromCell<T,DESCRIPTOR>(cell);
    return particle ? max((T)0.,(T)1.-particle->solidFraction) : (T)1.;
}

// 把孔隙率提取成一个标量场（后面算压降和孔隙体积都用它做权重）
class ExtractPorosity3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> blocks) {
        BlockLattice3D<T,DESCRIPTOR>& fluid = *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
        ScalarField3D<T>& eps = *dynamic_cast<ScalarField3D<T>*>(blocks[1]);
        const Dot3D offset = computeRelativeDisplacement(fluid,eps);
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z)
            eps.get(x+offset.x,y+offset.y,z+offset.z) = effectivePorosity(fluid.get(x,y,z));
    }
    virtual ExtractPorosity3D* clone() const { return new ExtractPorosity3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing; modified[1] = modif::staticVariables;
    }
};

// ── 二阶平衡态（全密度公式）。注意：这是普通函数，不能写 Descriptor<T>，要写 DESCRIPTOR<T> ──
static inline T scEquilibrium(plint iPop, T rho, Array<T,3> const& u) {
    T const cDotU = DESCRIPTOR<T>::c[iPop][0]*u[0] + DESCRIPTOR<T>::c[iPop][1]*u[1] + DESCRIPTOR<T>::c[iPop][2]*u[2];
    T const uSqr = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
    return DESCRIPTOR<T>::t[iPop]*rho*((T)1 + (T)3*cDotU + (T)4.5*cDotU*cDotU - (T)1.5*uSqr);
}
// 从格点读全密度和动量（Palabos 存的是 f_i − t_i，求和后 +1 才是密度）
static inline void cellRhoJ(Cell<T,DESCRIPTOR> const& cell, T& rho, Array<T,3>& j) {
    T rhoBar=0; j[0]=j[1]=j[2]=0;
    for (plint iPop=0; iPop<DESCRIPTOR<T>::q; ++iPop) {
        rhoBar += cell[iPop];
        for (int d=0; d<3; ++d) j[d] += cell[iPop]*DESCRIPTOR<T>::c[iPop][d];
    }
    rho = rhoBar + (T)1;
}

// =====================================================================================
//  入口条件：每步把入口面（x=0）的两套格子直接置为平衡态
//    mode=0 恒流量：速度给定 u_in，总密度取相邻格点（x=1）的值（零梯度，避免同时固定
//                  速度和压力造成的过约束——原版就把密度也固定成 1 了）
//    mode=1 恒压力：总密度给定 ρ_in（对应入口压力），速度取相邻格点的值
//    两种模式下入口组成都固定：注入相占 sInj，另一相 1−sInj
// =====================================================================================
class OilWaterInlet3D : public BoxProcessingFunctional3D_LL<T,DESCRIPTOR,T,DESCRIPTOR> {
public:
    OilWaterInlet3D(int mode_, T uIn_, T rhoIn_, T sInj_, bool waterInjected_)
        : mode(mode_), uIn(uIn_), rhoIn(rhoIn_), sInj(sInj_), waterInjected(waterInjected_) { }
    virtual void process(Box3D domain, BlockLattice3D<T,DESCRIPTOR>& latticeW, BlockLattice3D<T,DESCRIPTOR>& latticeO) {
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            // 相邻格点（下游一格）的总密度和总速度
            T rwN, roN; Array<T,3> jwN, joN;
            cellRhoJ(latticeW.get(x+1,y,z), rwN, jwN);
            cellRhoJ(latticeO.get(x+1,y,z), roN, joN);
            T const rhoTotN = max(rwN+roN, (T)1e-6);
            Array<T,3> u((T)0,(T)0,(T)0);
            T rhoTot;
            if (mode==0) { u[0] = uIn; rhoTot = rhoTotN; }          // 恒流量
            else { for (int d=0; d<3; ++d) u[d] = (jwN[d]+joN[d])/rhoTotN; rhoTot = rhoIn; }  // 恒压力
            T const rhoW_in = waterInjected ? sInj*rhoTot : ((T)1-sInj)*rhoTot;
            T const rhoO_in = waterInjected ? ((T)1-sInj)*rhoTot : sInj*rhoTot;
            Cell<T,DESCRIPTOR>& cellW = latticeW.get(x,y,z);
            Cell<T,DESCRIPTOR>& cellO = latticeO.get(x,y,z);
            for (plint iPop=0; iPop<DESCRIPTOR<T>::q; ++iPop) {
                cellW[iPop] = scEquilibrium(iPop,rhoW_in,u) - DESCRIPTOR<T>::t[iPop];   // 存成去偏移形式
                cellO[iPop] = scEquilibrium(iPop,rhoO_in,u) - DESCRIPTOR<T>::t[iPop];
            }
        }
    }
    virtual OilWaterInlet3D* clone() const { return new OilWaterInlet3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables; modified[1] = modif::staticVariables;
    }
private:
    int mode; T uIn, rhoIn, sInj; bool waterInjected;
};

// =====================================================================================
//  出口组成：出口面（x=nx−1）挂的是 Palabos 压力边界，它需要每个格点的目标密度。
//  原版把水、油同时设成本体密度（总密度 1.8），出口变成了双倍密度的全混合源（审计算例 2）。
//  这里：总密度固定为 rhoOut，水/油按上游一格的局部组成分配。
// =====================================================================================
class OilWaterOutletComposition3D : public BoxProcessingFunctional3D_LL<T,DESCRIPTOR,T,DESCRIPTOR> {
public:
    OilWaterOutletComposition3D(T rhoOut_, T rhoResidual_) : rhoOut(rhoOut_), rhoResidual(rhoResidual_) { }
    virtual void process(Box3D domain, BlockLattice3D<T,DESCRIPTOR>& latticeW, BlockLattice3D<T,DESCRIPTOR>& latticeO) {
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            T rwN, roN; Array<T,3> dummy;
            cellRhoJ(latticeW.get(x-1,y,z), rwN, dummy);
            cellRhoJ(latticeO.get(x-1,y,z), roN, dummy);
            T s = (rwN+roN>0) ? rwN/(rwN+roN) : (T)0.5;
            s = min(max(s, rhoResidual), (T)1-rhoResidual);
            latticeW.get(x,y,z).defineDensity(s*rhoOut);
            latticeO.get(x,y,z).defineDensity(((T)1-s)*rhoOut);
        }
    }
    virtual OilWaterOutletComposition3D* clone() const { return new OilWaterOutletComposition3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables; modified[1] = modif::staticVariables;
    }
private:
    T rhoOut, rhoResidual;
};

// =====================================================================================
//  两相标量场：含水饱和度 S_w = ρ_w/(ρ_w+ρ_o)（格子密度油水都为 1，直接按密度算）
//              两相总压 p_lb = c_s²(ρ_w+ρ_o) + c_s² G ρ_w ρ_o（Shan-Chen 状态方程）
// =====================================================================================
class OilWaterFields3D : public BoxProcessingFunctional3D {
public:
    OilWaterFields3D(T G_) : G(G_) { }
    virtual void processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> blocks) {
        BlockLattice3D<T,DESCRIPTOR>& latticeW = *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
        BlockLattice3D<T,DESCRIPTOR>& latticeO = *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[1]);
        ScalarField3D<T>& Sw = *dynamic_cast<ScalarField3D<T>*>(blocks[2]);
        ScalarField3D<T>& p  = *dynamic_cast<ScalarField3D<T>*>(blocks[3]);
        const Dot3D oS = computeRelativeDisplacement(latticeW,Sw);
        const Dot3D oP = computeRelativeDisplacement(latticeW,p);
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            T rw, ro; Array<T,3> dummy;
            cellRhoJ(latticeW.get(x,y,z), rw, dummy);
            cellRhoJ(latticeO.get(x,y,z), ro, dummy);
            T s = (rw+ro>0) ? rw/(rw+ro) : (T)0;
            Sw.get(x+oS.x,y+oS.y,z+oS.z) = min(max(s,(T)0),(T)1);
            p.get(x+oP.x,y+oP.y,z+oP.z) = ((T)1/(T)3)*(rw+ro) + ((T)1/(T)3)*G*rw*ro;
        }
    }
    virtual OilWaterFields3D* clone() const { return new OilWaterFields3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const {
        modified[0]=modif::nothing; modified[1]=modif::nothing;
        modified[2]=modif::staticVariables; modified[3]=modif::staticVariables;
    }
private:
    T G;
};

// =====================================================================================
//  ★★ 2026-10-25 新增：两相流场统计（只统计孔隙格点，即有效孔隙率 > 0.5 的格点）
//    · 总密度 ρ_w+ρ_o 的最大/最小值 —— Shan-Chen 的分相强度取决于 G×ρ_总，用来监测压差是否过大
//    · 混合物速度 |j_w+j_o|/(ρ_w+ρ_o) 的最大/平均值 —— 用来判断是否失稳
//    旧版用"水格子"自己的速度判失稳：油区里只有残余密度（0.005）的水，它的 j/ρ 没有物理意义、
//    在界面旁可以很大，会造成误报；现在改用两相合起来的速度。
//    这是 Palabos 的"归约型"处理器：每个进程统计自己那块，Palabos 自动做全局 max/sum。
// =====================================================================================
class TwoPhaseStats3D : public ReductiveBoxProcessingFunctional3D_LL<T,DESCRIPTOR,T,DESCRIPTOR> {
public:
    TwoPhaseStats3D()
        : maxRhoId(this->getStatistics().subscribeMax()), negMinRhoId(this->getStatistics().subscribeMax()),
          maxSpeedId(this->getStatistics().subscribeMax()), sumSpeedId(this->getStatistics().subscribeSum()),
          sumRhoWId(this->getStatistics().subscribeSum()), sumRhoOId(this->getStatistics().subscribeSum()),
          countId(this->getStatistics().subscribeSum()), nanId(this->getStatistics().subscribeSum()) { }
    virtual void process(Box3D domain, BlockLattice3D<T,DESCRIPTOR>& latticeW, BlockLattice3D<T,DESCRIPTOR>& latticeO) {
        BlockStatistics& st = this->getStatistics();
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            Cell<T,DESCRIPTOR>& cw = latticeW.get(x,y,z);
            if (effectivePorosity(cw) <= (T)0.5) continue;            // 壁面、颗粒内部不统计
            T rw, ro; Array<T,3> jw, jo;
            cellRhoJ(cw, rw, jw);
            cellRhoJ(latticeO.get(x,y,z), ro, jo);
            T const rt = rw + ro;
            T const ux = (jw[0]+jo[0])/rt, uy = (jw[1]+jo[1])/rt, uz = (jw[2]+jo[2])/rt;
            T const sp = std::sqrt(ux*ux + uy*uy + uz*uz);
            if (!std::isfinite((double)rt) || !std::isfinite((double)sp)) { st.gatherSum(nanId, 1.); continue; }
            st.gatherMax(maxRhoId, rt);   st.gatherMax(negMinRhoId, -rt);   // 最小值 = −max(−ρ)
            st.gatherMax(maxSpeedId, sp); st.gatherSum(sumSpeedId, sp);
            st.gatherSum(sumRhoWId, rw);  st.gatherSum(sumRhoOId, ro);  st.gatherSum(countId, 1.);
        }
    }
    virtual TwoPhaseStats3D* clone() const { return new TwoPhaseStats3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const { modified[0]=modif::nothing; modified[1]=modif::nothing; }
    T count() const       { return (T)this->getStatistics().getSum(countId); }
    T maxRho() const      { return (T)this->getStatistics().getMax(maxRhoId); }
    T minRho() const      { return -(T)this->getStatistics().getMax(negMinRhoId); }
    T maxSpeed() const    { return (T)this->getStatistics().getMax(maxSpeedId); }
    T meanSpeed() const   { T n=count(); return n>0 ? (T)this->getStatistics().getSum(sumSpeedId)/n : (T)0; }
    T meanRhoW() const    { T n=count(); return n>0 ? (T)this->getStatistics().getSum(sumRhoWId)/n : (T)0; }
    T meanRhoO() const    { T n=count(); return n>0 ? (T)this->getStatistics().getSum(sumRhoOId)/n : (T)0; }
    T nanCount() const    { return (T)this->getStatistics().getSum(nanId); }
private:
    plint maxRhoId, negMinRhoId, maxSpeedId, sumSpeedId, sumRhoWId, sumRhoOId, countId, nanId;
};

// =====================================================================================
//  颗粒周围含水饱和度采样（供胶结损伤律使用，P0-6 ①）：
//  遍历"部分被颗粒覆盖"的格点（0 < 固含率 ≤ SOLFRAC_MAX，即颗粒表面那一圈），
//  把该格的 ρ_w、ρ_o 累加到所属颗粒。结果数组按 LIGGGHTS 本地下标存放。
// =====================================================================================
class ParticleSaturationSampler3D : public BoxProcessingFunctional3D_LL<T,DESCRIPTOR,T,DESCRIPTOR> {
public:
    ParticleSaturationSampler3D(std::vector<T>* sumW_, std::vector<T>* sumO_, LiggghtsCouplingWrapper* wrapper_)
        : sumW(sumW_), sumO(sumO_), wrapper(wrapper_) { }
    virtual void process(Box3D domain, BlockLattice3D<T,DESCRIPTOR>& latticeW, BlockLattice3D<T,DESCRIPTOR>& latticeO) {
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            Cell<T,DESCRIPTOR>& cw = latticeW.get(x,y,z);
            IBdynamicsParticleData<T,DESCRIPTOR>* pd = getParticleDataFromCell<T,DESCRIPTOR>(cw);
            if (!pd || pd->partId<1 || pd->solidFraction<=SOLFRAC_MIN || pd->solidFraction>SOLFRAC_MAX) continue;
            plint ind = wrapper->lmp->atom->map((int)pd->partId);
            if (ind<0 || ind>=(plint)sumW->size()) continue;
            T rw, ro; Array<T,3> dummy;
            cellRhoJ(cw, rw, dummy);
            cellRhoJ(latticeO.get(x,y,z), ro, dummy);
            (*sumW)[ind] += max(rw,(T)0);
            (*sumO)[ind] += max(ro,(T)0);
        }
    }
    virtual ParticleSaturationSampler3D* clone() const { return new ParticleSaturationSampler3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& modified) const {
        modified[0]=modif::nothing; modified[1]=modif::nothing;
    }
private:
    std::vector<T>* sumW; std::vector<T>* sumO; LiggghtsCouplingWrapper* wrapper;
};

// ── VTK 输出：水/油速度、含水饱和度、两相总压（Pa，相对出口）────────────────────────
// =====================================================================================
//  ★★ 2026-10-26（提速）：固定颗粒（骨架 type1、筛网 type3）的映射结果只算一次、存起来。
//    LBDEMcoupling 原做法每步对"所有颗粒、两套格子"重新做一遍 setSpheresOnLattice：
//    每个骨架颗粒（直径约 30 格）表面附近的格点都要做 125 点子采样，N=20 测试中约占每步 25% 的时间，
//    而骨架和筛网根本不动，结果每步都一样。
//    现在：t=0 先只映射固定颗粒，把每个格点的（固含率、颗粒编号）存进两个标量场；
//    之后每步 ① 把两套格子的颗粒信息恢复成这份"静态底图"（颗粒速度为 0），② 只映射会动的细粒（type 2）。
//    细粒与骨架共享的格点仍按原规则"固含率大的颗粒占用"；细粒离开后该格点自动回到骨架的值
//    （原做法里这种格点会被清零，要靠骨架恰好排在后面重新映射才能恢复，结果依赖颗粒在数组里的顺序）。
// =====================================================================================
class StoreStaticParticles3D : public BoxProcessingFunctional3D {
public:
    // blocks：[0] 格子 [1] 固含率场 [2] 颗粒编号场
    virtual void processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> blocks) {
        BlockLattice3D<T,DESCRIPTOR>& lattice = *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
        ScalarField3D<T>& sf = *dynamic_cast<ScalarField3D<T>*>(blocks[1]);
        ScalarField3D<T>& id = *dynamic_cast<ScalarField3D<T>*>(blocks[2]);
        const Dot3D o1 = computeRelativeDisplacement(lattice, sf), o2 = computeRelativeDisplacement(lattice, id);
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            IBdynamicsParticleData<T,DESCRIPTOR>* pd = getParticleDataFromCell<T,DESCRIPTOR>(lattice.get(x,y,z));
            sf.get(x+o1.x,y+o1.y,z+o1.z) = pd ? pd->solidFraction : (T)0;
            id.get(x+o2.x,y+o2.y,z+o2.z) = pd ? (T)pd->partId : (T)0;
        }
    }
    virtual StoreStaticParticles3D* clone() const { return new StoreStaticParticles3D(*this); }
    virtual void getTypeOfModification(vector<modif::ModifT>& m) const { m[0]=modif::nothing; m[1]=modif::staticVariables; m[2]=modif::staticVariables; }
};

class RestoreStaticParticles3D : public BoxProcessingFunctional3D {
public:
    // blocks：[0] 水格子 [1] 油格子 [2] 固含率场 [3] 颗粒编号场
    virtual void processGenericBlocks(Box3D domain, vector<AtomicBlock3D*> blocks) {
        BlockLattice3D<T,DESCRIPTOR>& lw = *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[0]);
        BlockLattice3D<T,DESCRIPTOR>& lo = *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[1]);
        ScalarField3D<T>& sf = *dynamic_cast<ScalarField3D<T>*>(blocks[2]);
        ScalarField3D<T>& id = *dynamic_cast<ScalarField3D<T>*>(blocks[3]);
        const Dot3D oo = computeRelativeDisplacement(lw, lo);
        const Dot3D o1 = computeRelativeDisplacement(lw, sf), o2 = computeRelativeDisplacement(lw, id);
        const Array<T,3> zero((T)0,(T)0,(T)0);
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            const T s0 = sf.get(x+o1.x,y+o1.y,z+o1.z);
            const plint i0 = (plint)(id.get(x+o2.x,y+o2.y,z+o2.z)+(T)0.5);
            IBdynamicsParticleData<T,DESCRIPTOR>* pw = getParticleDataFromCell<T,DESCRIPTOR>(lw.get(x,y,z));
            if (pw) { pw->solidFraction = s0; pw->partId = i0; pw->uPart = zero; }
            IBdynamicsParticleData<T,DESCRIPTOR>* po = getParticleDataFromCell<T,DESCRIPTOR>(lo.get(x+oo.x,y+oo.y,z+oo.z));
            if (po) { po->solidFraction = s0; po->partId = i0; po->uPart = zero; }
        }
    }
    virtual RestoreStaticParticles3D* clone() const { return new RestoreStaticParticles3D(*this); }
    // 故意声明"不修改"（与 LBDEMcoupling 的 SetSingleSphere3D 同理）：紧接着的细粒映射最后会用
    // AttributeFunctional 统一同步一次各进程交界处的颗粒信息，这里不必多同步一次。
    virtual void getTypeOfModification(vector<modif::ModifT>& m) const { m[0]=modif::nothing; m[1]=modif::nothing; m[2]=modif::nothing; m[3]=modif::nothing; }
};

// 每步：颗粒 → 两套格子（恢复静态底图 + 只映射细粒）
static void setSpheresOnBothLattices(MultiBlockLattice3D<T,DESCRIPTOR>& latticeW, MultiBlockLattice3D<T,DESCRIPTOR>& latticeO,
                                     LiggghtsCouplingWrapper& wrapper, PhysUnits3D<T> const& units, Box3D const* couplingDomain,
                                     MultiScalarField3D<T>& staticSf, MultiScalarField3D<T>& staticId)
{
    vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&latticeO); a.push_back(&staticSf); a.push_back(&staticId);
    applyProcessingFunctional(new RestoreStaticParticles3D(), latticeW.getBoundingBox(), a);
    std::vector<plint> excludeFixed; excludeFixed.push_back(1); excludeFixed.push_back(3);   // 只映射 type 2 细粒
    setSpheresOnLattice(latticeW,wrapper,units,excludeFixed,false,couplingDomain);
    setSpheresOnLattice(latticeO,wrapper,units,excludeFixed,false,couplingDomain);
}

void writeTwoPhaseVTK(MultiBlockLattice3D<T,DESCRIPTOR>& latticeW, MultiBlockLattice3D<T,DESCRIPTOR>& latticeO,
                      PhysUnits3D<T> const& units, plint iter, T G, T pScale, T pRefLb, std::string const& prefix)
{
    std::string fname(createFileName(prefix+"vtk", iter, (int)environmentInteger("VTK_FILE_DIGITS")));
    VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
    vtkOut.writeData<3,float>(*computeVelocity(latticeW), "velocity_w", units.getPhysVel(1));
    vtkOut.writeData<3,float>(*computeVelocity(latticeO), "velocity_o", units.getPhysVel(1));
    MultiScalarField3D<T> Sw(latticeW), p(latticeW);
    vector<MultiBlock3D*> args; args.push_back(&latticeW); args.push_back(&latticeO); args.push_back(&Sw); args.push_back(&p);
    applyProcessingFunctional(new OilWaterFields3D(G), latticeW.getBoundingBox(), args);
    vtkOut.writeData<float>(Sw, "water_saturation", 1.);
    subtractInPlace(p, pRefLb);
    vtkOut.writeData<float>(p, "pressure_total", pScale);
    // ★ 追加到 .pvd 时间索引（ParaView 直接打开即可播动画）
    static bool pvdHeaderWritten = false;
    std::string pvdPath = prefix + "timeseries.pvd";
    {
        std::ofstream pvd(pvdPath, pvdHeaderWritten ? std::ios::app : std::ios::trunc);
        if (!pvdHeaderWritten) {
            pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            pvd << "  <Collection>\n";
            pvdHeaderWritten = true;
        }
        std::string base = fname;
        size_t slash = base.find_last_of("/\\");
        if (slash != std::string::npos) base = base.substr(slash+1);
        pvd << "    <DataSet timestep=\"" << iter << "\" group=\"\" part=\"0\" file=\""
            << base << ".vti\"/>\n";
    }
    pcout << "wrote " << fname << std::endl;
}

// ── 构造一套 PSM 复合动力学（base = 正则化 BGK + Smagorinsky（Cs=0 时退化为正则化 BGK））──
static IBcompositeDynamics<T,DESCRIPTOR>* makeIBDynamics(T omega, T cSmago) {
    return new IBcompositeDynamics<T,DESCRIPTOR>(new SecuredSmagorinskyRegularizedDynamics<T,DESCRIPTOR>(omega,cSmago));
}

// =====================================================================================
//  驱动协议：返回当前时刻的"驱动因子"（乘在入口速度或入口压差上）
//    continuous   : 0→1 线性缓启（rampTime），之后保持 1
//    stop_restart : 同上，直到注入达到 STOP_AT_PV；然后置 0 保持 HOLD_TIME_S（停注）；
//                   再分 RESTART_STAGES 级台阶升到 RESTART_MAX_FACTOR，每级保持 RESTART_STAGE_HOLD_S
//    stage 编号：−1 = 正常驱替，0 = 停注，1..N = 重启第几级
// =====================================================================================
struct DriveProtocol {
    bool stopRestart; T rampTime, holdTime, stageHold, maxFactor; plint stages; T stopAtPV;
    bool stopped; T tStop;     // 运行中状态：是否已进入停注，停注开始时刻
    DriveProtocol() : stopRestart(false), rampTime(1), holdTime(0), stageHold(0), maxFactor(1), stages(1), stopAtPV(1e30), stopped(false), tStop(0) {}
    T factor(T driveTime, T PV, int& stage) {
        stage = -1;
        if (!stopRestart) return min((T)1., driveTime/rampTime);
        if (!stopped && PV >= stopAtPV) { stopped = true; tStop = driveTime; }
        if (!stopped) return min((T)1., driveTime/rampTime);
        T const dt = driveTime - tStop;
        if (dt < holdTime) { stage = 0; return (T)0; }
        plint k = min<plint>((plint)std::floor((dt-holdTime)/stageHold) + 1, stages);
        stage = (int)k;
        return maxFactor*(T)k/(T)stages;
    }
};

// =====================================================================================
//                                            主函数
// =====================================================================================
int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    water_fangsha::ParameterConfig* parameterConfig = 0;
    try { parameterConfig = new water_fangsha::ParameterConfig(); }
    catch (std::exception const& error) { pcout << "[FATAL] 参数配置读取失败: " << error.what() << "\n"; return 1; }
    activeParameterConfig = parameterConfig;

    // ── 命令行参数（与单相/原版一致）──────────────────────────────────────────
    plint N(0); T uMax(0.), rho_f(0.), mu_f(0.), v_inf(0.), maxT(0.); std::string outDir;
    T nuScale = environmentDouble("NU_SCALE");
    try {
        global::argv(1).read(N); global::argv(2).read(uMax); global::argv(3).read(rho_f); global::argv(4).read(mu_f);
        global::argv(5).read(v_inf); global::argv(6).read(maxT); global::argv(7).read(outDir);
    } catch(PlbIOException&) {
        pcout << "用法: ./oilwaterDrive N uMax rho_f mu_f v_inf maxT outDir [nuScale dt_dem totalPV rampTime openOutletPV calibrationMode packFile kMinMd kMaxMd phiMin phiMax]\n";
        exit(1);
    }
    try { global::argv(8).read(nuScale); } catch(PlbIOException&) {}

    // ── 单相沿用的参数（含义见 defaults.env）──────────────────────────────────
    T rampTime = environmentDouble("RAMP_TIME");
    T dt_dem_target = environmentDouble("DT_DEM");
    try { global::argv(9).read(dt_dem_target); } catch(PlbIOException&) {}
    const T dSmall = environmentDouble("MIN_WEAK_DIAMETER_UM")*1.e-6;
    T openOutletPV = environmentDouble("OPEN_OUTLET_PV");
    T total_PV = environmentDouble("TOTAL_PV");
    try { global::argv(10).read(total_PV); } catch(PlbIOException&) {}
    try { global::argv(11).read(rampTime); } catch(PlbIOException&) {}
    try { global::argv(12).read(openOutletPV); } catch(PlbIOException&) {}
    plint calibrationMode = (plint)environmentInteger("CALIBRATION_MODE");
    try { global::argv(13).read(calibrationMode); } catch(PlbIOException&) {}
    std::string packFile = parameterConfig->text("PACK_FILE");
    std::string outputPrefix = parameterConfig->text("OUTPUT_PREFIX");
    try { global::argv(14).read(packFile); } catch(PlbIOException&) {}
    T initialKMinMd = environmentDouble("K_MIN_MD"), initialKMaxMd = environmentDouble("K_MAX_MD");
    T initialPhiMin = environmentDouble("PHI_MIN"), initialPhiMax = environmentDouble("PHI_MAX");
    try { global::argv(15).read(initialKMinMd); } catch(PlbIOException&) {}
    try { global::argv(16).read(initialKMaxMd); } catch(PlbIOException&) {}
    try { global::argv(17).read(initialPhiMin); } catch(PlbIOException&) {}
    try { global::argv(18).read(initialPhiMax); } catch(PlbIOException&) {}
    std::vector<T> pvSamples = environmentDoubleList("PV_SAMPLES");
    const T vtkT = environmentDouble("VTK_INTERVAL_S");
    const T particleDumpT = environmentDouble("PARTICLE_DUMP_INTERVAL_S");
    const T logT = environmentDouble("LOG_INTERVAL_S");
    const T smagorinskyConstant = environmentDouble("SMAGORINSKY_CONSTANT");
    const T steadyRelativeTolerance = environmentDouble("K_STEADY_REL_TOL");
    const T steadyWindowTime = environmentDouble("K_STEADY_WINDOW_S");
    const T steadyHoldTime = environmentDouble("K_STEADY_HOLD_S");
    const plint steadyWindowSize = max<plint>((plint)std::ceil(steadyWindowTime/logT)+1,3);
    const plint requiredSteadyCount = max<plint>((plint)std::ceil(steadyHoldTime/logT),1);
    const T minimumTimeSafetyFactor = environmentDouble("MIN_TIME_SAFETY_FACTOR");
    const T maxPredictedDensityDrop = environmentDouble("MAX_PREDICTED_LB_DENSITY_DROP");
    const T configuredEquilibrationTime = environmentDouble("EQUILIBRATION_TIME_S");
    const T demNeighborSkin = environmentDouble("DEM_NEIGHBOR_SKIN_M");
    const plint demNeighborDelay = environmentInteger("DEM_NEIGHBOR_DELAY_STEPS");
    const plint demNeighborEvery = environmentInteger("DEM_NEIGHBOR_EVERY_STEPS");
    const plint demNeighborOne = environmentInteger("DEM_NEIGHBOR_ONE");
    const plint demNeighborPage = environmentInteger("DEM_NEIGHBOR_PAGE");
    const T demYoungsModulus = environmentDouble("DEM_YOUNGS_MODULUS_PA");
    const T demPoissonRatio = environmentDouble("DEM_POISSON_RATIO");
    const T demRestitution = environmentDouble("DEM_RESTITUTION");
    const T demFriction = environmentDouble("DEM_FRICTION");
    const T screenYoungsModulus = environmentDouble("SCREEN_YOUNGS_MODULUS_PA");
    const T screenPoissonRatio = environmentDouble("SCREEN_POISSON_RATIO");
    const T screenRestitution = environmentDouble("SCREEN_RESTITUTION");
    const T screenFriction = environmentDouble("SCREEN_FRICTION");
    const T screenX = environmentDouble("SCREEN_X_M");
    const T screenCaptureBand = environmentDouble("SCREEN_CAPTURE_BAND_M");
    const T screenPressureOffset = environmentDouble("SCREEN_PRESSURE_OFFSET_M");
    const T screenBlockedSpeed = environmentDouble("SCREEN_BLOCKED_SPEED_M_S");
    T demMaxParticleSpeed = environmentDouble("DEM_MAX_PARTICLE_SPEED_M_S");   // 0 = 自动（按马赫数）
    const T demMaxParticleOmega = environmentDouble("DEM_MAX_PARTICLE_OMEGA_RAD_S");
    const T bondRadiusMultiplier = environmentDouble("BOND_RADIUS_MULTIPLIER");
    const T bondPfcEffectiveModulus = environmentDouble("BOND_PFC_EFFECTIVE_MODULUS_GPA")*1.0e9;
    const T bondPfcShearModulus = environmentDouble("BOND_PFC_SHEAR_MODULUS_GPA")*1.0e9;
    const T bondPfcStiffnessRatioInput = environmentDouble("BOND_PFC_STIFFNESS_RATIO");
    const T bondPfcStiffnessRatio = bondPfcStiffnessRatioInput>0. ? bondPfcStiffnessRatioInput : bondPfcEffectiveModulus/bondPfcShearModulus;
    const T bondSigma = environmentDouble("BOND_NORMAL_STRENGTH_MPA")*1.0e6;
    const T bondTau = environmentDouble("BOND_SHEAR_STRENGTH_MPA")*1.0e6;
    const T bondDnf = environmentDouble("BOND_DAMP_NORMAL_FORCE"), bondDtf = environmentDouble("BOND_DAMP_TANGENTIAL_FORCE");
    const T bondDnt = environmentDouble("BOND_DAMP_NORMAL_TORQUE"), bondDtt = environmentDouble("BOND_DAMP_TANGENTIAL_TORQUE");
    const T bondCreateDistance = environmentDouble("BOND_CREATE_DISTANCE_M");
    const T bondCreateRatio = environmentDouble("BOND_CREATE_RATIO");
    T ibmForceAbort = environmentDouble("IBM_FORCE_ABORT_N");
    T ibmTorqueAbort = environmentDouble("IBM_TORQUE_ABORT_NM");
    const T demViscousDamping = environmentDouble("DEM_VISCOUS_DAMPING");
    const plint demThermoEvery = environmentInteger("DEM_THERMO_EVERY");
    const plint particleCleanupInterval = environmentInteger("PARTICLE_CLEANUP_INTERVAL_LB_STEPS");
    const T minParticleResolutionCells = environmentDouble("MIN_PARTICLE_RESOLUTION_CELLS");
    const T minLbVelocity = environmentDouble("MIN_LB_VELOCITY");
    const T maxLbVelocity = environmentDouble("MAX_LB_VELOCITY");
    const plint envelopeWidth = environmentInteger("LBM_ENVELOPE_WIDTH");
    const plint ibOutletExcludedLayers = environmentInteger("IBM_OUTLET_EXCLUDED_LAYERS");
    const plint diagnosticScanSteps = environmentInteger("DIAGNOSTIC_SCAN_STEPS");   // ★★ 前 N 步扫描异常格点（0=关）
    const plint pressurePlaneThickness = environmentInteger("PRESSURE_PLANE_THICKNESS_CELLS");
    const plint screenPressureHalfWidth = environmentInteger("SCREEN_PRESSURE_HALF_WIDTH_CELLS");

    // ── 两相参数（本版新增/修改的在注释里标 ★）────────────────────────────────
    const T rho_o = environmentDouble("RHO_O");              // 油密度 [kg/m³]，只用于油的运动粘度和 Bond 数
    const T mu_o  = environmentDouble("MU_O");               // 油动力粘度 [Pa·s]；粘度比 M = mu_o/mu_f
    const T shanChenG = environmentDouble("SHAN_CHEN_G");    // 跨组分排斥强度（验证：G=2.0~3.0 稳定，3.5 发散）
    const T rhoResidual = environmentDouble("RHO_RESIDUAL"); // 另一相残余密度（格子单位）
    std::string initialContinuousPhase = environmentText("INITIAL_CONTINUOUS_PHASE");
    std::string injectPhase = environmentText("INJECT_PHASE");
    const T initialFrontPosition = environmentDouble("INITIAL_FRONT_POSITION_M");
    const T tauWTarget = environmentDouble("TAU_W_TARGET");                    // ★ 水相松弛时间（≥0.8）
    // ★★ 2026-10-26：润湿性改用"几何接触角"模型，参数直接给接触角（度，从水相一侧量：<90 水湿，>90 油湿）
    Array<T,3> wallContactAngle(environmentDouble("WALL_CONTACT_ANGLE_DEG_TYPE1"),   // ★ 骨架
                                environmentDouble("WALL_CONTACT_ANGLE_DEG_TYPE2"),   // ★ 细粒
                                environmentDouble("WALL_CONTACT_ANGLE_DEG_TYPE3"));  // ★ 筛网
    const T containerWallContactAngle = environmentDouble("CONTAINER_WALL_CONTACT_ANGLE_DEG");   // 域侧壁（默认 90 中性）
    const T wettingGain = environmentDouble("WETTING_GAIN");                       // 接触角外推系数（标定值）
    const T frontWidthCells = environmentDouble("INITIAL_FRONT_WIDTH_CELLS");   // ★ 初始界面 tanh 半宽（格）
    const plint capillaryOnParticles = environmentInteger("CAPILLARY_FORCE_ON_PARTICLES"); // ★ 1=传毛细反作用力
    std::string inletControl = environmentText("INLET_CONTROL");                // ★ flow | pressure
    const T inletPressurePa = environmentDouble("INLET_PRESSURE_PA");           // ★ 恒压模式的入口压差
    std::string driveProtocolName = environmentText("DRIVE_PROTOCOL");          // ★ continuous | stop_restart
    const T stopAtPV = environmentDouble("STOP_AT_PV");
    const T holdTimeS = environmentDouble("HOLD_TIME_S");
    const plint restartStages = environmentInteger("RESTART_STAGES");
    const T restartStageHoldS = environmentDouble("RESTART_STAGE_HOLD_S");
    const T restartMaxFactor = environmentDouble("RESTART_MAX_FACTOR");
    const T particleMachLimit = environmentDouble("PARTICLE_MACH_LIMIT");       // ★ 细粒速度上限（格子马赫数）
    const T sigmaLbCalibrated = environmentDouble("SIGMA_LB_CALIBRATED");       // ★ 由 validation/v2_laplace.py 标定
    const T sigmaPhysReference = environmentDouble("SIGMA_PHYS_REFERENCE_N_M"); // 真实油水界面张力（仅用于提示）
    // ★★ 2026-10-25 新增（修正超算报错时加入，详见 修改说明_20261025.md）
    const plint mixtureTau = environmentInteger("MIXTURE_TAU");          // 1 = 按局部含水比例混合粘度（粘度比 M 较大时必须开）
    const T mixSatLow  = environmentDouble("MIX_SAT_LOW");               // 局部含水比例低于它 → 纯油粘度
    const T mixSatHigh = environmentDouble("MIX_SAT_HIGH");              // 局部含水比例高于它 → 纯水粘度
    const T traceVelocityCap = environmentDouble("TRACE_VELOCITY_CAP");  // 痕量组分"偏离公共速度"的上限（格子单位）
    const T diffusionTau = environmentDouble("DIFFUSION_TAU");          // ★★ 2026-10-26：组分间动量交换松弛时间 τ_D（0=旧做法）
    const T scSolidThreshold = environmentDouble("SHAN_CHEN_SOLID_THRESHOLD");   // ★★ 2026-10-26：固含率超过它的格点在 Shan-Chen 里当固体
    const T maxRuntimeDensityRise = environmentDouble("MAX_RUNTIME_LB_DENSITY_RISE"); // 运行中实测入口−出口密度差超过它就停算

    bool pvSamplesValid = true;
    for (size_t i=0; i<pvSamples.size(); ++i) if (pvSamples[i]<=0. || (i>0 && pvSamples[i]<=pvSamples[i-1])) pvSamplesValid=false;
    if (pvSamplesValid && calibrationMode==0) {
        bool targetIncluded=false; std::vector<T>::iterator insertion = pvSamples.begin();
        for (;insertion!=pvSamples.end() && *insertion<total_PV;++insertion) if(fabs(*insertion-total_PV)<(T)1.e-12) targetIncluded=true;
        if(insertion!=pvSamples.end() && fabs(*insertion-total_PV)<(T)1.e-12) targetIncluded=true;
        if(!targetIncluded) pvSamples.insert(insertion,total_PV);
    }
    const bool pressureControl = (inletControl=="pressure");
    const bool stopRestart = (driveProtocolName=="stop_restart");

    // ── 参数合法性检查（2026-10-25 重写）──────────────────────────────────────────
    //   旧版把 70 多个条件写在同一个 if 里，任何一个不满足都只打印同一句
    //   "参数无效…（oilwater_fangsha 要求 OPEN_OUTLET_PV=0）"，
    //   让人误以为是 OPEN_OUTLET_PV 的问题（你报错记录里的"错误 4"就是这样）。
    //   现在逐条检查：哪个参数不合法、当前值是多少、要求是什么，全部打印出来再退出。
    int nParamErrors = 0;
    auto requireParam = [&](bool ok, const char* name, double value, const char* rule) {
        if (!ok) { ++nParamErrors; pcout << "[PARAM_ERROR] " << name << " = " << value << "  → 要求：" << rule << "\n"; }
    };
    auto requireText = [&](bool ok, const char* name, std::string const& value, const char* rule) {
        if (!ok) { ++nParamErrors; pcout << "[PARAM_ERROR] " << name << " = '" << value << "'  → 要求：" << rule << "\n"; }
    };
    // 网格、流体、时间
    requireParam(N>0, "N", (double)N, "> 0（每个参考直径 REFERENCE_DIAMETER_UM 分成几个格子）");
    requireParam(uMax>0., "UMAX", uMax, "> 0（TAU_W_TARGET>0 时会被自动改写，但仍要求为正）");
    requireParam(rho_f>0., "RHO_F", rho_f, "> 0");
    requireParam(mu_f>0., "MU_F", mu_f, "> 0");
    requireParam(v_inf>0., "V_INF", v_inf, "> 0");
    requireParam(maxT>0., "MAX_T", maxT, "> 0");
    requireParam(nuScale>0., "NU_SCALE", nuScale, "> 0");
    requireParam(dt_dem_target>0., "DT_DEM", dt_dem_target, "> 0");
    requireParam(rampTime>0., "RAMP_TIME", rampTime, "> 0");
    requireParam(openOutletPV==(T)0., "OPEN_OUTLET_PV", openOutletPV, "必须 = 0（本算例从一开始就回收出口细粒）");
    requireParam(total_PV>0., "TOTAL_PV", total_PV, "> 0");
    requireParam(calibrationMode==0 || calibrationMode==1, "CALIBRATION_MODE", (double)calibrationMode, "0（生产）或 1（单相渗透率标定）");
    requireParam(initialKMinMd>0., "K_MIN_MD", initialKMinMd, "> 0");
    requireParam(initialKMaxMd>initialKMinMd, "K_MAX_MD", initialKMaxMd, "> K_MIN_MD");
    requireParam(initialPhiMin>0., "PHI_MIN", initialPhiMin, "> 0");
    requireParam(initialPhiMax>initialPhiMin && initialPhiMax<1., "PHI_MAX", initialPhiMax, "PHI_MIN < PHI_MAX < 1");
    requireParam(dSmall>0., "MIN_WEAK_DIAMETER_UM", dSmall*1e6, "> 0");
    requireParam(vtkT>0., "VTK_INTERVAL_S", vtkT, "> 0");
    requireParam(particleDumpT>0., "PARTICLE_DUMP_INTERVAL_S", particleDumpT, "> 0");
    requireParam(logT>0., "LOG_INTERVAL_S", logT, "> 0");
    requireParam(smagorinskyConstant>=0., "SMAGORINSKY_CONSTANT", smagorinskyConstant, ">= 0（建议 0）");
    requireParam(steadyRelativeTolerance>0. && steadyRelativeTolerance<1., "K_STEADY_REL_TOL", steadyRelativeTolerance, "0 < 值 < 1");
    requireParam(steadyWindowTime>0., "K_STEADY_WINDOW_S", steadyWindowTime, "> 0");
    requireParam(steadyHoldTime>0., "K_STEADY_HOLD_S", steadyHoldTime, "> 0");
    requireParam(minimumTimeSafetyFactor>=(T)1., "MIN_TIME_SAFETY_FACTOR", minimumTimeSafetyFactor, ">= 1");
    requireParam(maxPredictedDensityDrop>0. && maxPredictedDensityDrop<=(T)0.35, "MAX_PREDICTED_LB_DENSITY_DROP", maxPredictedDensityDrop,
                 "0 < 值 <= 0.35（推荐 0.1–0.2）。它是入口与出口之间允许的格子总密度差；"
                 "Shan-Chen 的分相强度取决于 G×ρ_总，ρ_总 变化太大界面张力会沿岩心漂移，G×ρ_总>3.4 会发散。"
                 "想加快注入请看 修改说明_20261025.md 第 4 问，不要放宽这个值");
    requireParam(maxRuntimeDensityRise>0. && maxRuntimeDensityRise<=(T)0.5, "MAX_RUNTIME_LB_DENSITY_RISE", maxRuntimeDensityRise,
                 "0 < 值 <= 0.5（推荐 0.35）：运行中实测的入口−出口密度差超过它就停算，避免把算力花在失真的结果上");
    requireParam(configuredEquilibrationTime>=0., "EQUILIBRATION_TIME_S", configuredEquilibrationTime, ">= 0");
    // DEM 邻居表与接触参数
    requireParam(demNeighborSkin>0., "DEM_NEIGHBOR_SKIN_M", demNeighborSkin, "> 0");
    requireParam(demNeighborDelay>=0, "DEM_NEIGHBOR_DELAY_STEPS", (double)demNeighborDelay, ">= 0");
    requireParam(demNeighborEvery>0, "DEM_NEIGHBOR_EVERY_STEPS", (double)demNeighborEvery, "> 0");
    requireParam(demNeighborOne>0, "DEM_NEIGHBOR_ONE", (double)demNeighborOne, "> 0");
    requireParam(demNeighborPage>0, "DEM_NEIGHBOR_PAGE", (double)demNeighborPage, "> 0");
    requireParam(demYoungsModulus>0., "DEM_YOUNGS_MODULUS_PA", demYoungsModulus, "> 0");
    requireParam(demPoissonRatio>=0. && demPoissonRatio<0.5, "DEM_POISSON_RATIO", demPoissonRatio, "0 <= 值 < 0.5");
    requireParam(demRestitution>=0. && demRestitution<=1., "DEM_RESTITUTION", demRestitution, "0 <= 值 <= 1");
    requireParam(demFriction>=0., "DEM_FRICTION", demFriction, ">= 0");
    requireParam(screenYoungsModulus>0., "SCREEN_YOUNGS_MODULUS_PA", screenYoungsModulus, "> 0");
    requireParam(screenPoissonRatio>=0. && screenPoissonRatio<0.5, "SCREEN_POISSON_RATIO", screenPoissonRatio, "0 <= 值 < 0.5");
    requireParam(screenRestitution>=0. && screenRestitution<=1., "SCREEN_RESTITUTION", screenRestitution, "0 <= 值 <= 1");
    requireParam(screenFriction>=0., "SCREEN_FRICTION", screenFriction, ">= 0");
    requireParam(screenX>0., "SCREEN_X_M", screenX, "> 0");
    requireParam(screenCaptureBand>0., "SCREEN_CAPTURE_BAND_M", screenCaptureBand, "> 0");
    requireParam(screenPressureOffset>0., "SCREEN_PRESSURE_OFFSET_M", screenPressureOffset, "> 0");
    requireParam(screenBlockedSpeed>=0., "SCREEN_BLOCKED_SPEED_M_S", screenBlockedSpeed, ">= 0");
    requireParam(demMaxParticleSpeed>=0., "DEM_MAX_PARTICLE_SPEED_M_S", demMaxParticleSpeed, ">= 0（0 = 按格子马赫数自动）");
    requireParam(demMaxParticleOmega>0., "DEM_MAX_PARTICLE_OMEGA_RAD_S", demMaxParticleOmega, "> 0");
    // 胶结键
    requireParam(bondRadiusMultiplier>0., "BOND_RADIUS_MULTIPLIER", bondRadiusMultiplier, "> 0");
    requireParam(bondPfcEffectiveModulus>0., "BOND_PFC_EFFECTIVE_MODULUS_GPA", bondPfcEffectiveModulus/1e9, "> 0");
    requireParam(bondSigma>0., "BOND_NORMAL_STRENGTH_MPA", bondSigma/1e6, "> 0");
    requireParam(bondTau>0., "BOND_SHEAR_STRENGTH_MPA", bondTau/1e6, "> 0");
    requireParam(bondCreateDistance>0., "BOND_CREATE_DISTANCE_M", bondCreateDistance, "> 0");
    requireParam(bondCreateRatio>0., "BOND_CREATE_RATIO", bondCreateRatio, "> 0");
    requireParam(demViscousDamping>=0., "DEM_VISCOUS_DAMPING", demViscousDamping, ">= 0");
    requireParam(demThermoEvery>0, "DEM_THERMO_EVERY", (double)demThermoEvery, "> 0");
    requireParam(particleCleanupInterval>0, "PARTICLE_CLEANUP_INTERVAL_LB_STEPS", (double)particleCleanupInterval, "> 0");
    requireParam(pvSamplesValid, "PV_SAMPLES", pvSamples.empty() ? -1. : (double)pvSamples[0], "逗号分隔、全部 > 0、严格递增，例如 1,5,10");
    // 数值检查阈值与统计区域
    requireParam(minParticleResolutionCells>0., "MIN_PARTICLE_RESOLUTION_CELLS", minParticleResolutionCells, "> 0");
    requireParam(minLbVelocity>0., "MIN_LB_VELOCITY", minLbVelocity, "> 0");
    requireParam(maxLbVelocity>minLbVelocity, "MAX_LB_VELOCITY", maxLbVelocity, "> MIN_LB_VELOCITY");
    requireParam(envelopeWidth>0, "LBM_ENVELOPE_WIDTH", (double)envelopeWidth, "> 0");
    requireParam(ibOutletExcludedLayers>=1, "IBM_OUTLET_EXCLUDED_LAYERS", (double)ibOutletExcludedLayers,
                 ">= 1（两相相互作用不能算到出口面上，否则会读到出口外未初始化的格点；推荐 2）");
    requireParam(pressurePlaneThickness>0, "PRESSURE_PLANE_THICKNESS_CELLS", (double)pressurePlaneThickness, "> 0");
    requireParam(screenPressureHalfWidth>=0, "SCREEN_PRESSURE_HALF_WIDTH_CELLS", (double)screenPressureHalfWidth, ">= 0");
    // 两相
    requireParam(rho_o>0., "RHO_O", rho_o, "> 0");
    requireParam(mu_o>0., "MU_O", mu_o, "> 0");
    requireParam(shanChenG>=0., "SHAN_CHEN_G", shanChenG, ">= 0（验证过 2.0–3.0；3.5 发散）");
    requireParam(rhoResidual>0. && rhoResidual<0.5, "RHO_RESIDUAL", rhoResidual, "0 < 值 < 0.5（推荐 0.005）");
    requireText(initialContinuousPhase=="oil" || initialContinuousPhase=="water", "INITIAL_CONTINUOUS_PHASE", initialContinuousPhase, "oil 或 water");
    requireText(injectPhase=="oil" || injectPhase=="water", "INJECT_PHASE", injectPhase, "oil 或 water");
    requireParam(initialFrontPosition>=0., "INITIAL_FRONT_POSITION_M", initialFrontPosition, ">= 0");
    requireParam(tauWTarget>0.5, "TAU_W_TARGET", tauWTarget,
                 "> 0.5（推荐 0.7–0.8；numpy 验证中 0.6 的小液滴/接触角算例发散）。"
                 "两相版本不再支持 0（沿用 UMAX 会得到 τ≈0.518，Shan-Chen 必然发散）");
    requireParam(wallContactAngle[0]>=20. && wallContactAngle[0]<=160., "WALL_CONTACT_ANGLE_DEG_TYPE1", wallContactAngle[0], "20 <= 值 <= 160（度）");
    requireParam(wallContactAngle[1]>=20. && wallContactAngle[1]<=160., "WALL_CONTACT_ANGLE_DEG_TYPE2", wallContactAngle[1], "20 <= 值 <= 160（度）");
    requireParam(wallContactAngle[2]>=20. && wallContactAngle[2]<=160., "WALL_CONTACT_ANGLE_DEG_TYPE3", wallContactAngle[2], "20 <= 值 <= 160（度）");
    requireParam(containerWallContactAngle>=20. && containerWallContactAngle<=160., "CONTAINER_WALL_CONTACT_ANGLE_DEG", containerWallContactAngle, "20 <= 值 <= 160（度，推荐 90）");
    requireParam(wettingGain>0. && wettingGain<=5., "WETTING_GAIN", wettingGain, "0 < 值 <= 5（标定值，一般不用改）");
    // 润湿系数 K = gain·cotθ（θ=90° 时为 0：壁面对流体不施任何附加力）
    Array<T,3> wallK;
    for (int k=0;k<3;++k) { const T r = wallContactAngle[k]*(T)3.141592653589793/(T)180; wallK[k] = wettingGain*std::cos(r)/std::sin(r); }
    const T containerWallK = wettingGain*std::cos(containerWallContactAngle*(T)3.141592653589793/(T)180)
                                        /std::sin(containerWallContactAngle*(T)3.141592653589793/(T)180);
    requireParam(frontWidthCells>=0., "INITIAL_FRONT_WIDTH_CELLS", frontWidthCells, ">= 0（推荐 4）");
    requireText(inletControl=="flow" || inletControl=="pressure", "INLET_CONTROL", inletControl, "flow 或 pressure");
    requireParam(!pressureControl || inletPressurePa>0., "INLET_PRESSURE_PA", inletPressurePa, "恒压模式下 > 0");
    requireText(driveProtocolName=="continuous" || driveProtocolName=="stop_restart", "DRIVE_PROTOCOL", driveProtocolName, "continuous 或 stop_restart");
    if (stopRestart) {
        requireParam(stopAtPV>0., "STOP_AT_PV", stopAtPV, "> 0（stop_restart 协议）");
        requireParam(holdTimeS>=0., "HOLD_TIME_S", holdTimeS, ">= 0（stop_restart 协议）");
        requireParam(restartStages>0, "RESTART_STAGES", (double)restartStages, "> 0（stop_restart 协议）");
        requireParam(restartStageHoldS>0., "RESTART_STAGE_HOLD_S", restartStageHoldS, "> 0（stop_restart 协议）");
        requireParam(restartMaxFactor>0., "RESTART_MAX_FACTOR", restartMaxFactor, "> 0（stop_restart 协议）");
    }
    requireParam(particleMachLimit>0. && particleMachLimit<=(T)0.2, "PARTICLE_MACH_LIMIT", particleMachLimit, "0 < 值 <= 0.2（推荐 0.05；PSM 要求颗粒每步移动远小于 1 格）");
    requireParam(sigmaLbCalibrated>=0., "SIGMA_LB_CALIBRATED", sigmaLbCalibrated, ">= 0（用 validation/v2_laplace.py 标定）");
    requireParam(mixtureTau==0 || mixtureTau==1, "MIXTURE_TAU", (double)mixtureTau, "0 或 1（推荐 1）");
    requireParam(mixSatLow>=0. && mixSatLow<mixSatHigh && mixSatHigh<=1., "MIX_SAT_LOW/MIX_SAT_HIGH", mixSatLow,
                 "0 <= MIX_SAT_LOW < MIX_SAT_HIGH <= 1（推荐 0.05 / 0.95）");
    requireParam(traceVelocityCap>0. && traceVelocityCap<=(T)0.5, "TRACE_VELOCITY_CAP", traceVelocityCap, "0 < 值 <= 0.5（推荐 0.3）");
    requireParam(diffusionTau==0. || (diffusionTau>(T)0.5 && diffusionTau<=(T)2.), "DIFFUSION_TAU", diffusionTau, "0（旧做法）或 0.5 < 值 <= 2（推荐 1.0）");
    requireParam(scSolidThreshold>=(T)0.3 && scSolidThreshold<=(T)0.999, "SHAN_CHEN_SOLID_THRESHOLD", scSolidThreshold, "0.3 <= 值 <= 0.999（推荐 0.5）");
    if (nParamErrors>0) {
        pcout << "[FATAL] 共 " << nParamErrors << " 个参数不合法（见上面的 [PARAM_ERROR] 行）。"
              << "请在 config/defaults.env（或 CONFIG_FILE 指定的剖面文件）里修改后重新提交。\n";
        return 1;
    }

    // ── 输出目录 ─────────────────────────────────────────────────────────
    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    if (global::mpi().getRank() == 0) {
        std::string mkdirCmd = "mkdir -p \"" + outDir + "\" \"" + lbOutDir + "\" \"" + demOutDir + "\"";
        if (std::system(mkdirCmd.c_str()) != 0) pcout << "[WARN] mkdir 返回非零\n";
    }
    global::mpi().barrier();
    global::directories().setOutputDir(lbOutDir);

    // ── 物理参数与单位换算 ───────────────────────────────────────────────
    const T nu_phys = mu_f/rho_f;              // 水运动粘度
    const T nu_f = nuScale*nu_phys;
    const T refD = environmentDouble("REFERENCE_DIAMETER_UM")*1.e-6;
    const T lx = environmentDouble("DOMAIN_LX_M"), ly = environmentDouble("DOMAIN_LY_M"), lz = environmentDouble("DOMAIN_LZ_M");
    const T matrixXLo = environmentDouble("CORE_X0_M"), matrixXHi = environmentDouble("CORE_X1_M");
    const T demBoxXHi = environmentDouble("DEM_BOX_XHI_M");
    if (!(refD>0. && lx>0. && ly>0. && lz>0. && matrixXLo>=0. && matrixXHi>matrixXLo && matrixXHi<screenX && screenX<lx
          && screenPressureOffset<screenX-matrixXLo && screenPressureOffset<lx-screenX && demBoxXHi>=lx)) {
        pcout << "[FATAL] 域尺寸/核心区/参考直径环境参数无效\n"; return 1;
    }
    const T r_ = refD/2.0;

    // ★ P0-1：若给定 TAU_W_TARGET，则反算 uMax 使水相 τ 等于目标值。
    //   Palabos 的 IncomprFlowParam 里 ν_lb = uMax·N/Re，τ = 0.5 + 3ν_lb，所以
    //   uMax = (τ−0.5)/3 · Re / N，其中 Re = refD·v_inf/ν_f。
    //   原版沿用单相的 UMAX=0.0015 得到 τ≈0.518，Shan-Chen 在这个 τ 下不可能稳定。
    const T Re_phys = refD*v_inf/nu_f;
    if (tauWTarget > 0.) {
        const T uMaxAuto = (tauWTarget-(T)0.5)/(T)3*Re_phys/(T)N;
        pcout << "[TAU] TAU_W_TARGET=" << tauWTarget << " → uMax 由 " << uMax << " 改为 " << uMaxAuto << "\n";
        uMax = uMaxAuto;
    }

    PhysUnits3D<T> units(2.*r_, v_inf, nu_f, lx, ly, lz, N, uMax, rho_f);
    units.setLbOffset(0.,0.,0.);
    IncomprFlowParam<T> parameters(units.getLbParam());
    const T physDx = units.getPhysLength(1);
    const T dt_phys = units.getPhysTime(1);

    // ── 把 DEM 需要的变量交给 LIGGGHTS 脚本 ───────────────────────────────
    LiggghtsCouplingWrapper wrapper(argv, MPI_COMM_WORLD);
    wrapper.setVariable("rho_fluid", rho_f);
    wrapper.setVariable("dmp_dir", demOutDir);
    wrapper.setVariable("pack_file", packFile);
    wrapper.setVariable("output_prefix", outputPrefix);
    wrapper.setVariable("dem_neighbor_skin",demNeighborSkin);
    wrapper.setVariable("dem_neighbor_delay",(T)demNeighborDelay);
    wrapper.setVariable("dem_neighbor_every",(T)demNeighborEvery);
    wrapper.setVariable("dem_neighbor_one",(T)demNeighborOne);
    wrapper.setVariable("dem_neighbor_page",(T)demNeighborPage);
    wrapper.setVariable("dem_youngs_modulus",demYoungsModulus);
    wrapper.setVariable("dem_poisson_ratio",demPoissonRatio);
    wrapper.setVariable("dem_restitution",demRestitution);
    wrapper.setVariable("dem_friction",demFriction);
    wrapper.setVariable("screen_youngs_modulus",screenYoungsModulus);
    wrapper.setVariable("screen_poisson_ratio",screenPoissonRatio);
    wrapper.setVariable("screen_restitution",screenRestitution);
    wrapper.setVariable("screen_friction",screenFriction);
    // ★ P0-7：细粒速度上限。0 表示自动：格子马赫数 PARTICLE_MACH_LIMIT × (Δx/Δt)。
    //   原版固定 0.05 m/s，比平均孔隙流速（0.05/0.365≈0.14 m/s）还低，再运移被人为截断。
    if (demMaxParticleSpeed <= 0.) demMaxParticleSpeed = particleMachLimit*physDx/dt_phys;
    // ★★ 2026-10-26：手动给的上限若超过 0.2 格/步，PSM 会失真（本机测试中颗粒 0.3 格/步时周围流场发散）
    if (demMaxParticleSpeed > (T)0.2*physDx/dt_phys)
        pcout << "[WARN] DEM_MAX_PARTICLE_SPEED_M_S=" << demMaxParticleSpeed << " m/s 超过 0.2 格/步（" << (T)0.2*physDx/dt_phys
              << " m/s）：颗粒移动太快时 PSM 失真、可能发散，建议设 0 用 PARTICLE_MACH_LIMIT=0.05 自动计算。\n";
    pcout << "[DEM] 细粒速度上限 = " << demMaxParticleSpeed << " m/s（" << demMaxParticleSpeed*dt_phys/physDx << " 格/步）\n";
    wrapper.setVariable("dem_max_particle_speed",demMaxParticleSpeed);
    wrapper.setVariable("dem_max_particle_omega",demMaxParticleOmega);
    wrapper.setVariable("bond_radius_multiplier",bondRadiusMultiplier);
    wrapper.setVariable("bond_pfc_effective_modulus",bondPfcEffectiveModulus);
    wrapper.setVariable("bond_pfc_stiffness_ratio",bondPfcStiffnessRatio);
    wrapper.setVariable("bond_sigma",bondSigma);
    wrapper.setVariable("bond_tau",bondTau);
    wrapper.setVariable("bond_dnf",bondDnf); wrapper.setVariable("bond_dtf",bondDtf);
    wrapper.setVariable("bond_dnt",bondDnt); wrapper.setVariable("bond_dtt",bondDtt);
    wrapper.setVariable("bond_create_distance",bondCreateDistance);
    wrapper.setVariable("bond_create_ratio",bondCreateRatio);
    wrapper.setVariable("dem_viscous_damping",demViscousDamping);
    wrapper.setVariable("dem_thermo_every",(T)demThermoEvery);
    wrapper.setVariable("phys_dx", physDx);
    wrapper.setVariable("domain_lx", lx); wrapper.setVariable("domain_ly", ly); wrapper.setVariable("domain_lz", lz);
    wrapper.setVariable("dem_box_xhi", demBoxXHi);
    wrapper.setVariable("outlet_x", lx);

    // ── 第一个 LIGGGHTS 脚本：建初始颗粒体系（几何与相无关）─────────────────
    wrapper.execFile("in.lbdem");

    double solidLocal=0., weakCountLocal0=0., skeletonCountLocal0=0., screenCountLocal0=0., weakRadiusMinLocal=1e9;
    for (int ia=0; ia<wrapper.lmp->atom->nlocal; ++ia) {
        double rr = wrapper.lmp->atom->radius[ia];
        solidLocal += 4.0*3.141592653589793*rr*rr*rr/3.0;
        if (wrapper.lmp->atom->type[ia]==2) { weakCountLocal0 += 1.; weakRadiusMinLocal = min(weakRadiusMinLocal, rr); }
        if (wrapper.lmp->atom->type[ia]==1) skeletonCountLocal0 += 1.;
        if (wrapper.lmp->atom->type[ia]==3) screenCountLocal0 += 1.;
    }
    double solid0=0., weakCount0=0., skeletonCount0=0., screenCount0=0., weakRadiusMin=1e9;
    MPI_Allreduce(&solidLocal,&solid0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&weakCountLocal0,&weakCount0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&skeletonCountLocal0,&skeletonCount0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&screenCountLocal0,&screenCount0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&weakRadiusMinLocal,&weakRadiusMin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    const T particleVolumePorosity = max((T)0.,min((T)1.,(T)(1.0-solid0/(lx*ly*lz))));
    pcout << "[GEOMETRY_DIAGNOSTIC] wholeBoxParticleVolumePorosity=" << particleVolumePorosity
          << " initialSkeletonParticles=" << skeletonCount0 << " initialWeakParticles=" << weakCount0
          << " fixedScreenParticles=" << screenCount0 << std::endl;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint particleDumpSteps = max<plint>(units.getLbSteps(particleDumpT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "oilwater_fangsha_youhua: oil-water two-phase LBM-DEM fines migration");
    pcout << "grid : " << parameters.getNx() << " × " << parameters.getNy() << " × " << parameters.getNz() << std::endl;

    LatticeDecomposition lDec(parameters.getNx(), parameters.getNy(), parameters.getNz(), wrapper.lmp);
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();

    // ── 两个弛豫频率：水从参数；油按运动粘度比换算（ν_o/ν_w = (μ_o/μ_w)(ρ_w/ρ_o)）──
    const T omegaW = parameters.getOmega();
    const T nu_w_lb = ((T)1/omegaW - (T)0.5)/(T)3;
    const T nu_o_lb = nu_w_lb*((mu_o/rho_o)/(mu_f/rho_f));
    const T tau_o = (T)0.5 + (T)3*nu_o_lb;
    const T omegaO = (T)1/tau_o;
    // ★ P0-2：两相格子本体密度都取 1。Shan-Chen 用理想气体状态方程，格子密度差会变成
    //   伪压力（原版 ρ_o=0.8 相当于约 74 kPa 的压力跳变）。真实油水密度差的影响用 Bond 数衡量。
    const T rhoTotalBulk = (T)1.0;   // 本体总密度（格子单位）：水区、油区都取 1

    // ── 无量纲数（写进方法节用）────────────────────────────────────────────
    const T phiGuess = 0.5*(initialPhiMin+initialPhiMax);
    const T poreVelocity = v_inf/phiGuess;
    const T M_visc = mu_o/mu_f;
    const T Bo = fabs(rho_f-rho_o)*9.81*(dSmall/2)*(dSmall/2)/sigmaPhysReference;
    const T pScale = rho_f*(physDx/dt_phys)*(physDx/dt_phys);      // 1 格子压力单位 = 多少 Pa
    const T sigmaScale = rho_f*physDx*physDx*physDx/(dt_phys*dt_phys); // 1 格子张力单位 = 多少 N/m
    const T sigmaPhysModel = sigmaLbCalibrated*sigmaScale;
    const T Ca_model = (sigmaLbCalibrated>0.) ? mu_f*poreVelocity/sigmaPhysModel : -1.;
    const T Ca_real = mu_f*poreVelocity/sigmaPhysReference;
    pcout << "═══════ 两相参数与无量纲数 ═══════\n"
          << "  Δx=" << physDx*1e6 << " μm  Δt=" << dt_phys << " s  uIn_lb=" << v_inf/units.getPhysVel(1) << "\n"
          << "  τ_w=" << 1/omegaW << " τ_o=" << tau_o << "  粘度比 M=" << M_visc << "  Re_pore≈" << poreVelocity*dSmall/nu_phys << "\n"
          << "  G=" << shanChenG << "  σ_lb(标定)=" << sigmaLbCalibrated << " → σ_model=" << sigmaPhysModel*1e3 << " mN/m（真实≈" << sigmaPhysReference*1e3 << "）\n"
          << "  Ca_model=" << Ca_model << "  Ca_real=" << Ca_real << "  Bo=" << Bo << "\n"
          << "  接触角(type1,2,3)=" << wallContactAngle[0] << "," << wallContactAngle[1] << "," << wallContactAngle[2] << "° 侧壁=" << containerWallContactAngle << "°\n"
          << "  入口控制=" << inletControl << "  协议=" << driveProtocolName << "  细粒速度上限=" << demMaxParticleSpeed << " m/s\n"
          << "═══════════════════════════════════\n";
    // ★★ 2026-10-25：τ_o 的检查按"是否用局部混合松弛时间"分开说明（旧版只有一句笼统的警告）
    if (mixtureTau==1) {
        pcout << "[TAU] MIXTURE_TAU=1：每个格点按局部含水比例混合粘度（纯水区 τ=" << 1/omegaW << "，纯油区 τ=" << tau_o
              << "，含水比例在 " << mixSatLow << "–" << mixSatHigh << " 之间调和平均）\n";
        if (tau_o > 10.) pcout << "[WARN] τ_o=" << tau_o << " 很大：油区松弛很慢，正则化 BGK 精度下降"
                               << "（numpy 分层流验证到 τ_o=19 仍稳定；τ_o≈8 时总流量误差约 3%，但粘的那层流量偏大约 35%），"
                               << "而且每注 1 PV 的步数与粘度比成正比。\n";
    } else {
        if (tau_o >= 4.5) {
            pcout << "[FATAL] MIXTURE_TAU=0（两套格子各用固定 τ）且 τ_o=" << tau_o << " ≥ 4.5：有壁面的验证里 20 步内发散。\n"
                  << "        请设 MIXTURE_TAU=1（推荐），或降低 MU_O / TAU_W_TARGET 使 τ_o<3.5。\n";
            return 1;
        }
        if (1/omegaW < (T)0.75) {
            pcout << "[FATAL] MIXTURE_TAU=0（两套格子各用固定 τ）且 τ_w=" << 1/omegaW << " < 0.75：numpy 分层流验证中 τ_w=0.7 时"
                  << "即使粘度比只有 2.5 也发散（旧做法只在 τ_w≈0.8 验证过）。请设 MIXTURE_TAU=1，或把 TAU_W_TARGET 改为 0.8。\n";
            return 1;
        }
        if (tau_o > 3.5) pcout << "[WARN] MIXTURE_TAU=0 且 τ_o=" << tau_o << " > 3.5：接近发散边界，且水区溶入的油会让水的有效粘度明显偏大；建议 MIXTURE_TAU=1。\n";
    }
    if (Ca_model > 0. && fabs(log10(Ca_model/Ca_real)) > 0.5)
        pcout << "[WARN] 模型毛细数与真实毛细数相差超过 3 倍：σ、τ、Δx 绑定在一起，若要匹配 Ca 请调 V_INF 或 G。\n";
    // 胶结键量级检查（P0-6 ②）：键法向破坏力 vs 细粒的 Stokes 拖曳与毛细力量级
    {
        const T Rb = bondRadiusMultiplier*weakRadiusMin;
        const T Fbreak = bondSigma*3.141592653589793*Rb*Rb;
        const T Fdrag = 3*3.141592653589793*mu_f*dSmall*poreVelocity;
        const T Fcap = 2*3.141592653589793*(dSmall/2)*sigmaPhysReference;
        pcout << "[BOND_SCALE] 键破坏力≈" << Fbreak << " N，Stokes 拖曳≈" << Fdrag << " N，毛细力量级≈" << Fcap << " N\n";
        if (Fbreak > 100*max(Fdrag,Fcap)) pcout << "[WARN] 键强比流体/毛细载荷高两个数量级以上：流动不可能冲断胶结，请重新标定 BOND_NORMAL_STRENGTH_MPA。\n";
    }

    // ── 建两套格子（同一块分解）────────────────────────────────────────────
    MultiBlockLattice3D<T,DESCRIPTOR> latticeW(MultiBlockManagement3D(blockStructure,threadAttribution,envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(), defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(), makeIBDynamics(omegaW,smagorinskyConstant));
    defineDynamics(latticeW, latticeW.getBoundingBox(), makeIBDynamics(omegaW,smagorinskyConstant));
    MultiBlockLattice3D<T,DESCRIPTOR> latticeO(MultiBlockManagement3D(lDec.getBlockDistribution(),lDec.getThreadAttribution(),envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(), defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(), makeIBDynamics(omegaO,smagorinskyConstant));
    defineDynamics(latticeO, latticeO.getBoundingBox(), makeIBDynamics(omegaO,smagorinskyConstant));
    for (int d=0; d<3; ++d) { latticeW.periodicity().toggle(d,false); latticeO.periodicity().toggle(d,false); }

    const plint demSubsteps = max<plint>((plint)std::ceil(dt_phys/dt_dem_target),1);
    const T resSmall = dSmall/physDx;
    pcout << "[CHECK] 最小弱颗粒分辨率=" << resSmall << " cells/d（审稿建议 ≥20，当前默认 10，请做收敛性对比）\n";
    if (resSmall < minParticleResolutionCells) { pcout << "[FATAL] 最小弱颗粒不足 " << minParticleResolutionCells << " 格\n"; return 1; }
    const T uInLb = v_inf/units.getPhysVel(1);
    if (uInLb < minLbVelocity) pcout << "[WARN] uInLb=" << uInLb << " 过小\n";
    if (uInLb > maxLbVelocity) pcout << "[WARN] uInLb=" << uInLb << " 过大（可压缩误差）\n";

    // ── 边界条件 ───────────────────────────────────────────────────────────
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* bcW = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* bcO = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    const plint nx=parameters.getNx(), ny=parameters.getNy(), nz=parameters.getNz();
    Box3D inlet(0,0,1,ny-2,1,nz-2);
    Box3D outlet(nx-1,nx-1,1,ny-2,1,nz-2);
    Box3D yLo(0,nx-1,0,0,1,nz-2), yHi(0,nx-1,ny-1,ny-1,1,nz-2), zLo(0,nx-1,0,ny-1,0,0), zHi(0,nx-1,0,ny-1,nz-1,nz-1);
    // 入口：自定义平衡态处理器（每步覆盖）；出口：局部正则化压力边界；四侧壁：反弹
    bcW->setPressureConditionOnBlockBoundaries(latticeW, outlet);
    bcO->setPressureConditionOnBlockBoundaries(latticeO, outlet);
    // （四个侧壁的 BounceBack 改到下面"初始相分布"之后再定义，原因见那里的注释）
    const Array<T,3> uZero((T)0.,(T)0.,(T)0.);

    // ── 初始相分布（★ P0-3：tanh 平滑界面，半宽 frontWidthCells 格）────────────
    //   标定模式：全域水。生产模式：x<front 为注入相一侧（入口段），x≥front 为初始连续相（岩心）。
    const T frontLb = (calibrationMode==1) ? (T)(nx+10) : initialFrontPosition/physDx;
    T rwOutlet0 = (T)0.5*rhoTotalBulk, roOutlet0 = (T)0.5*rhoTotalBulk;   // 出口面的初始组成（循环里记下）
    for (plint x=0; x<nx; ++x) {
        T phiCore;   // 本切片中"岩心初始连续相"的体积分数 0..1
        if (calibrationMode==1) phiCore = 0.;
        else if (frontWidthCells>0.) phiCore = 0.5*(1.+std::tanh(((T)x-frontLb)/frontWidthCells));
        else phiCore = ((T)x>=frontLb) ? 1. : 0.;
        T const sw = (initialContinuousPhase=="oil") ? 1.-phiCore : phiCore;          // 含水分数
        // ★★ 2026-10-26：界面过渡区按"Shan-Chen 压力处处相等"给总密度。
        //   旧版过渡区两相各 0.5、总密度 1，可 Shan-Chen 压力 p = c_s²(ρ_w+ρ_o) + c_s²Gρ_wρ_o 在那里
        //   比本体高 0.2（格子单位，约 150 Pa），第 0 步就向两边各发出一个压力脉冲（云端复现：入口段
        //   回流速度 0.09、贴近界面的颗粒受力冲高）。真实 Shan-Chen 界面里总密度本来就是凹下去的，
        //   这里直接解 ρ_t + G s(1−s) ρ_t² = P0（P0 = 本体的同一组合）得到过渡区总密度 ρ_t，
        //   再按含水比例 s 分给两相：本体仍是 ρ_总=1、残余相 RHO_RESIDUAL。
        T const s = rhoResidual + ((T)1-(T)2*rhoResidual)*sw;                           // 含水比例（残余相保底）
        T const aSC = shanChenG*s*((T)1-s);
        T const P0 = rhoTotalBulk + shanChenG*rhoResidual*((T)1-rhoResidual)*rhoTotalBulk*rhoTotalBulk;
        T const rhoT = (aSC>(T)1e-12) ? (std::sqrt((T)1+(T)4*aSC*P0)-(T)1)/((T)2*aSC) : P0;
        T const rw = s*rhoT;
        T const ro = ((T)1-s)*rhoT;
        initializeAtEquilibrium(latticeW, Box3D(x,x,0,ny-1,0,nz-1), rw, uZero);
        initializeAtEquilibrium(latticeO, Box3D(x,x,0,ny-1,0,nz-1), ro, uZero);
        if (x==nx-1) { rwOutlet0 = rw; roOutlet0 = ro; }
    }
    // ★★ 2026-10-26 修订：侧壁的 BounceBack 必须在"按初始组成置平衡态"之后才定义。
    //   旧版先定义 BounceBack 再 initializeAtEquilibrium：BounceBack 的"平衡态"是 0（Palabos 存 f−t），
    //   侧壁格点里于是装着"油、水密度都为 1、静止"的分布函数。第 0 步反弹时这些分布函数被原样
    //   送进贴壁流体：油区贴壁格点凭空多出约 0.17 的水（水区多出油），并获得约 0.14 格子单位的
    //   离壁速度——云端复现中这正是第 0 步贴壁格点速度 0.15、总密度 1.2–1.3 的来源。
    //   现在侧壁格点先和旁边流体一样按本切片组成置成静止平衡态，再换成 BounceBack，
    //   第 0 步反弹出来的就是"和流体同组成、静止"的分布函数。
    defineDynamics(latticeW, yLo, new BounceBack<T,DESCRIPTOR>((T)1.)); defineDynamics(latticeW, yHi, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(latticeW, zLo, new BounceBack<T,DESCRIPTOR>((T)1.)); defineDynamics(latticeW, zHi, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(latticeO, yLo, new BounceBack<T,DESCRIPTOR>((T)1.)); defineDynamics(latticeO, yHi, new BounceBack<T,DESCRIPTOR>((T)1.));
    defineDynamics(latticeO, zLo, new BounceBack<T,DESCRIPTOR>((T)1.)); defineDynamics(latticeO, zHi, new BounceBack<T,DESCRIPTOR>((T)1.));
    // 出口压力边界的初始目标密度 = 出口面的初始组成（旧版两相各取 0.5，第 0 帧 VTK 出口面显示成"油水各半"）
    setBoundaryDensity(latticeW, outlet, rwOutlet0);
    setBoundaryDensity(latticeO, outlet, roOutlet0);
    latticeW.initialize(); latticeO.initialize();

    MultiScalarField3D<T> porosityField(latticeW);
    setToConstant(porosityField,porosityField.getBoundingBox(),(T)0.);
    // ★★ 2026-10-26：Shan-Chen 内核每步先把两相密度和固体标记抽到这三个工作场里，读邻居时直接查表（提速）
    MultiScalarField3D<T> rhoWField(latticeW), rhoOField(latticeW), solidFlagField(latticeW);
    const T dt_dem = dt_phys/(T)demSubsteps;
    wrapper.setVariable("t_step", dt_dem);
    wrapper.setVariable("dmp_stp", particleDumpSteps*demSubsteps);
    wrapper.execFile("in2.lbdem");
    const int integratorFixIndex = wrapper.lmp->modify->find_fix("integr_weak");
    if (integratorFixIndex<0) { pcout << "[FATAL] 找不到 integr_weak fix\n"; return 1; }
    LAMMPS_NS::Fix* integratorFix = wrapper.lmp->modify->fix[integratorFixIndex];
    wrapper.runUpto(demSubsteps-1);
    const int bondComputeIndex = wrapper.lmp->modify->find_compute("bond_status");
    if (bondComputeIndex<0) { pcout << "[FATAL] 找不到 bond_status compute\n"; return 1; }
    LAMMPS_NS::Compute* bondCompute = wrapper.lmp->modify->compute[bondComputeIndex];
    bondCompute->compute_vector();
    double initialActiveBonds = global::mpi().getRank()==0 ? bondCompute->vector[2] : 0.;
    MPI_Bcast(&initialActiveBonds,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    pcout << "[BOND_AUDIT] initialActiveBonds=" << initialActiveBonds << std::endl;


    // ── 各类统计区域 ───────────────────────────────────────────────────────
    Box3D ibCouplingDomain(0, nx-1-ibOutletExcludedLayers, 0, ny-1, 0, nz-1);
    Box3D shanChenDomain(1, nx-1-ibOutletExcludedLayers, 1, ny-2, 1, nz-2);   // 流体内部
    Box3D fluxPlane(1,1,1,ny-2,1,nz-2);
    const plint matrixX0 = max<plint>((plint)floor(matrixXLo/physDx+(T)0.5),1);
    const plint matrixX1 = min<plint>((plint)floor(matrixXHi/physDx+(T)0.5),nx-2);
    Box3D matrixCore(matrixX0,matrixX1,1,ny-2,1,nz-2);
    Box3D pressureLeft(matrixX0,min(matrixX0+pressurePlaneThickness-1,matrixX1),1,ny-2,1,nz-2);
    Box3D pressureRight(max(matrixX1-pressurePlaneThickness+1,matrixX0),matrixX1,1,ny-2,1,nz-2);
    const plint screenUpX = max<plint>((plint)floor((screenX-screenPressureOffset)/physDx+(T)0.5),matrixX1+1);
    const plint screenDownX = min<plint>((plint)floor((screenX+screenPressureOffset)/physDx+(T)0.5),nx-2);
    Box3D screenPressureUp(max<plint>(screenUpX-screenPressureHalfWidth,matrixX1+1),min<plint>(screenUpX+screenPressureHalfWidth,screenDownX-1),1,ny-2,1,nz-2);
    Box3D screenPressureDown(max<plint>(screenDownX-screenPressureHalfWidth,screenUpX+1),min<plint>(screenDownX+screenPressureHalfWidth,nx-2),1,ny-2,1,nz-2);
    Box3D totalPressureIn(1,min<plint>(2,nx-2),1,ny-2,1,nz-2);
    Box3D totalPressureOut(max<plint>(nx-3,1),nx-2,1,ny-2,1,nz-2);

    // t=0 映射固相、算孔隙率
    //   ★★ 2026-10-26：先只映射固定颗粒（骨架、筛网），存成"静态底图"，再按每步的做法叠加细粒（见 setSpheresOnBothLattices）
    MultiScalarField3D<T> staticSolidFraction(latticeW), staticParticleId(latticeW);
    { std::vector<plint> excludeFines(1,(plint)2);
      setSpheresOnLattice(latticeW,wrapper,units,excludeFines,false,&ibCouplingDomain);
      vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&staticSolidFraction); a.push_back(&staticParticleId);
      applyProcessingFunctional(new StoreStaticParticles3D(), latticeW.getBoundingBox(), a); }
    setSpheresOnBothLattices(latticeW,latticeO,wrapper,units,&ibCouplingDomain,staticSolidFraction,staticParticleId);
    { vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&porosityField);
      applyProcessingFunctional(new ExtractPorosity3D(),latticeW.getBoundingBox(),a); }
    const T channelCrossSectionNodes = (T)(ny-2)*(T)(nz-2);
    const T matrixGeometricNodes = channelCrossSectionNodes*(T)(matrixX1-matrixX0+1);
    const T initialPoreVolumeLb = computeSum(porosityField,matrixCore);
    const T phi0 = initialPoreVolumeLb/matrixGeometricNodes;
    if (!(phi0>=initialPhiMin && phi0<=initialPhiMax)) {
        pcout << (calibrationMode==0 ? "[FATAL_GEOMETRY]" : "[CALIBRATION_GEOMETRY_WARNING]") << " squareCorePorosity=" << phi0
              << " outside accepted range " << initialPhiMin << ".." << initialPhiMax << std::endl;
        pcout << "        说明：这是颗粒包在当前网格（N=" << N << "，Δx=" << physDx*1e6 << " μm）上离散后的岩心孔隙度，"
              << "会随 N 略变；若只是略低于 PHI_MIN（如 0.3476 vs 0.35），把 PHI_MIN 调到略低于它即可。\n";
        if (calibrationMode==0) return 1;
    }
    const T channelAreaPhys = channelCrossSectionNodes*physDx*physDx;
    const T initialPoreVolumePhys = initialPoreVolumeLb*physDx*physDx*physDx;
    const T pressureLengthPhys = (T)max<plint>(matrixX1-matrixX0-2,1)*physDx;
    const T pressurePerDensity = pScale/(T)3.;     // 1 个格子密度单位对应多少 Pa：ρ_f·c_s²·(Δx/Δt)²

    // ── 压差预检（PRESSURE_PREFLIGHT，2026-10-25 重写）──────────────────────────────
    //   LBM 里压力就是密度：p = c_s²ρ。岩心两端压差 Δp 必须靠入口比出口"更密"来实现，
    //       Δρ = Δp / (ρ_f · c_s² · (Δx/Δt)²)
    //   恒流量时用达西定律估 Δp = V_INF · μ · L / k_min（L 为岩心长度，k_min 为渗透率下限）。
    //   ★ 旧版的 bug：μ 用的是水的粘度，但生产算例的岩心一开始充满油（INITIAL_CONTINUOUS_PHASE=oil），
    //     开始驱替时的真实压差要再大 μ_o/μ_w 倍（M=20 时 20 倍）。现在取油、水粘度中较大者。
    //   Δρ 不能太大：Shan-Chen 的分相强度取决于 G×ρ_总，出口 ρ_总=1、入口 ρ_总=1+Δρ，
    //   Δρ 太大时界面张力沿岩心变化、G×ρ_总>3.4 还会发散（v1 验证 G=3.5 发散）。
    const T muWorst = (calibrationMode==1) ? mu_f : max(mu_f, mu_o);
    const T kMinM2 = initialKMinMd*(T)9.869233e-16;
    // stop_restart 协议最后一级会把驱动放大到 RESTART_MAX_FACTOR 倍，压差也跟着放大
    const T driveFactorMax = (stopRestart && calibrationMode==0) ? max((T)1., restartMaxFactor) : (T)1.;
    const T dpWaterFilled = v_inf*mu_f*pressureLengthPhys/kMinM2;
    const T predictedWorstDp = driveFactorMax*(pressureControl ? inletPressurePa : v_inf*muWorst*pressureLengthPhys/kMinM2);
    const T predictedDensityDrop = predictedWorstDp/pressurePerDensity;
    // 估算"多快才合格"以及"要算多少步"：注 1 PV 的物理时间 = 孔隙体积 / (达西速度 × 截面积)
    //   恒流：允许的最大流速；恒压：允许的最大入口压差 pressureMaxAllowed 对应的（最粘相）达西流速
    const T pressureMaxAllowed = maxPredictedDensityDrop*pressurePerDensity/driveFactorMax;
    const T vMaxAllowed = pressureControl ? pressureMaxAllowed*kMinM2/(muWorst*pressureLengthPhys)
                                          : v_inf*maxPredictedDensityDrop/predictedDensityDrop;
    const T vNow = pressureControl ? inletPressurePa*kMinM2/(muWorst*pressureLengthPhys) : v_inf;   // 恒压：按最粘的相估流速
    const T timePerPV = initialPoreVolumePhys/(vNow*channelAreaPhys);
    const T timePerPVAtVMax = initialPoreVolumePhys/(vMaxAllowed*channelAreaPhys);
    const T cellsTotal = (T)nx*(T)ny*(T)nz;
    pcout << "[PRESSURE_PREFLIGHT] 1 个格子密度单位 = " << pressurePerDensity << " Pa（Δx=" << physDx*1e6 << " μm, Δt=" << dt_phys << " s）\n"
          << "[PRESSURE_PREFLIGHT] 岩心充满水时 Δp≈" << dpWaterFilled << " Pa → Δρ≈" << dpWaterFilled/pressurePerDensity
          << "；充满" << (mu_o>mu_f ? "油" : "水") << "时（最不利）Δp≈" << predictedWorstDp << " Pa → Δρ≈" << predictedDensityDrop
          << (driveFactorMax>1. ? "（已乘重启最大倍数 RESTART_MAX_FACTOR）" : "")
          << "；上限 MAX_PREDICTED_LB_DENSITY_DROP=" << maxPredictedDensityDrop << std::endl;
    pcout << "[PRESSURE_PREFLIGHT] predicted_dp_Pa=" << predictedWorstDp << " predicted_delta_rho=" << predictedDensityDrop << " limit=" << maxPredictedDensityDrop << std::endl;
    if (calibrationMode==0 && predictedDensityDrop>maxPredictedDensityDrop) {
        const T stepsPerPV = timePerPVAtVMax/dt_phys;
        const T needT = minimumTimeSafetyFactor*(configuredEquilibrationTime + total_PV*timePerPVAtVMax + (T)0.5*rampTime);
        pcout << "[FATAL] 预测的入口−出口格子密度差 Δρ=" << predictedDensityDrop << " 超过上限 " << maxPredictedDensityDrop << "。\n"
              << "        在当前 N、TAU_W_TARGET、MU_O 下，" << (pressureControl ? "入口压差必须 INLET_PRESSURE_PA ≤ " : "流速必须 V_INF ≤ ")
              << (pressureControl ? pressureMaxAllowed : vMaxAllowed) << (pressureControl ? " Pa" : " m/s") << "；\n"
              << "        那时每注 1 PV 约需 " << stepsPerPV << " 个 LBM 步，注 TOTAL_PV=" << total_PV << " PV 需 MAX_T ≥ " << needT
              << " s（≈" << needT/dt_phys << " 步，每步 " << cellsTotal << " 个格点×2 套格子）。\n"
              << "        注意：加大 N 不能减少每 PV 的步数（见 修改说明_20261025.md 第 3 问）；能减少步数的是降低 TAU_W_TARGET、\n"
              << "        降低粘度比（MU_O）或缩短岩心。提交前可先运行  python3 plan_run.py  估算。\n";
        return 1;
    }
    const T equilibrationTime = calibrationMode==1 ? (T)0. : configuredEquilibrationTime;
    T idealDriveTime=0., requiredMaxTime=0.;
    if (calibrationMode==0 && !pressureControl) {
        const T idealFlowRate = v_inf*channelAreaPhys, targetVolume = total_PV*initialPoreVolumePhys, rampVolume = (T)0.5*idealFlowRate*rampTime;
        idealDriveTime = targetVolume<=rampVolume ? sqrt((T)2.*targetVolume*rampTime/idealFlowRate) : targetVolume/idealFlowRate+(T)0.5*rampTime;
        if (stopRestart) idealDriveTime += holdTimeS + restartStages*restartStageHoldS;
        requiredMaxTime = minimumTimeSafetyFactor*(equilibrationTime+idealDriveTime);
        pcout << "[TIME_BUDGET] 注 1 PV 物理时间≈" << timePerPV << " s（" << timePerPV/dt_phys << " 步）；TOTAL_PV=" << total_PV
              << " 需要 MAX_T ≥ " << requiredMaxTime << " s ≈ " << requiredMaxTime/dt_phys << " 步；当前 MAX_T=" << maxT
              << " s（" << maxT/dt_phys << " 步）\n";
        if (maxT<requiredMaxTime) {
            pcout << "[FATAL] MAX_T=" << maxT << " s 不够：以 V_INF=" << v_inf << " m/s 注入 TOTAL_PV=" << total_PV << " PV，"
                  << "需要 " << idealDriveTime << " s 驱替 + " << equilibrationTime << " s 平衡，再乘安全系数 " << minimumTimeSafetyFactor << "。\n"
                  << "        请把 MAX_T 设为 ≥ " << requiredMaxTime << "（或减小 TOTAL_PV）。这一步不改变物理结果，只是给够计算时长。\n";
            return 1;
        }
    } else if (pressureControl) {
        pcout << "[INFO] 恒压模式下流量未知，不做 maxT 预检；以 MAX_T 或 TOTAL_PV 先到者为准。按最粘的相估计：注 1 PV 约 "
              << timePerPV << " s（" << timePerPV/dt_phys << " 步）。\n";
    }
    if (global::mpi().getRank()==0) {
        ofstream geometryCsv((outDir+"geometry_metrics.csv").c_str());
        geometryCsv << "core_porosity,accepted_porosity_min,accepted_porosity_max,core_x0_m,core_x1_m,cross_section_y_m,cross_section_z_m,channel_area_m2,initial_pore_volume_m3,whole_box_particle_volume_porosity,screen_x_m,initial_screen_particles,screen_pressure_up_x_m,screen_pressure_down_x_m\n";
        geometryCsv << setprecision(17) << phi0 << "," << initialPhiMin << "," << initialPhiMax << "," << matrixX0*physDx << "," << matrixX1*physDx << ","
                    << (ny-2)*physDx << "," << (nz-2)*physDx << "," << channelAreaPhys << "," << initialPoreVolumePhys << "," << particleVolumePorosity << ","
                    << screenX << "," << screenCount0 << "," << screenUpX*physDx << "," << screenDownX*physDx << "\n";
        ofstream timeCsv((outDir+"time_budget.csv").c_str());
        timeCsv << "calibration_mode,target_PV,equilibration_s,ideal_drive_s,safety_factor,required_max_time_s,configured_max_time_s\n";
        timeCsv << setprecision(17) << calibrationMode << "," << total_PV << "," << equilibrationTime << "," << idealDriveTime << "," << minimumTimeSafetyFactor << "," << requiredMaxTime << "," << maxT << "\n";
        ofstream dimless((outDir+"dimensionless_numbers.csv").c_str());
        dimless << "dx_m,dt_s,tau_w,tau_o,viscosity_ratio,G,sigma_lb,sigma_model_N_m,sigma_real_N_m,Ca_model,Ca_real,Bo,Re_pore,particle_speed_limit_m_s\n";
        dimless << setprecision(10) << physDx << "," << dt_phys << "," << 1/omegaW << "," << tau_o << "," << M_visc << "," << shanChenG << "," << sigmaLbCalibrated << ","
                << sigmaPhysModel << "," << sigmaPhysReference << "," << Ca_model << "," << Ca_real << "," << Bo << "," << poreVelocity*dSmall/nu_phys << "," << demMaxParticleSpeed << "\n";
    }

    // ── 入口组成与压力 ───────────────────────────────────────────────────
    const bool waterInjected = (injectPhase=="water");
    const T sInj = (T)1 - rhoResidual;                              // 注入相在入口的质量分数
    const T inletDeltaRho = inletPressurePa/pressurePerDensity;     // 恒压模式：入口总密度 = 出口 + Δρ·factor
    DriveProtocol protocol;
    protocol.stopRestart = stopRestart && calibrationMode==0; protocol.rampTime = rampTime; protocol.holdTime = holdTimeS;
    protocol.stageHold = restartStageHoldS; protocol.maxFactor = restartMaxFactor; protocol.stages = restartStages; protocol.stopAtPV = stopAtPV;

    // ★ 两相力中转站 + 颗粒饱和度采样数组
    TwoPhaseForceBuffers<T> capBuffers;
    std::vector<T> satSumW, satSumO;
    std::map<int,double> immersionTime;     // 颗粒编号 → 累计处于水相（S_w>0.5）的时间 [s]

    // ===================== 主时间推进循环 =====================
    T injectedVolume=0., PV=0.; size_t nextPV=0;
    ofstream permCsv, calibrationCsv, forceCsv, balanceCsv, bondCsv, flowCsv, screenCsv, protocolCsv, satCsv;
    if (global::mpi().getRank()==0) {
        permCsv.open((outDir+"perm_porosity.csv").c_str());
        permCsv << "PV,target_PV,t_s,core_porosity,total_mobility_m2_Pa_s,apparent_permeability_md_water_visc,core_dp_Pa,total_dp_Pa,q_water_m3_s,q_oil_m3_s,n_produced,n_remaining\n";
        screenCsv.open((outDir+"screen_metrics.csv").c_str());
        screenCsv << "lb_step,t_s,PV,screen_dp_Pa,total_dp_Pa,near_screen_upstream,near_screen_slow,downstream_of_screen,produced_weak\n";
        if (calibrationMode==1) {
            calibrationCsv.open((outDir+"permeability_calibration.csv").c_str());
            calibrationCsv << "lb_step,t_s,core_porosity,flow_rate_m3_s,core_dp_Pa,total_dp_Pa,permeability_m2,permeability_md,instant_relative_change,window_relative_range,steady_count\n";
        } else {
            forceCsv.open((outDir+"ibm_force_audit.csv").c_str());
            forceCsv << "lb_step,t_s,max_force_N,max_torque_Nm,max_capillary_force_N,max_particle_speed_m_s,max_particle_surface_speed_m_s,velocity_limit_count,over_force,over_torque,nonfinite_count,max_fixed_particle_force_N\n";
            balanceCsv.open((outDir+"particle_balance.csv").c_str());
            balanceCsv << "lb_step,t_s,initial_skeleton,remaining_skeleton,initial_screen,remaining_screen,initial_weak,remaining_weak,produced_weak,abnormal_loss\n";
            bondCsv.open((outDir+"bond_history.csv").c_str());
            bondCsv << "lb_step,t_s,initial_active_bonds,active_bonds,broken_this_lb_step,cumulative_broken_bonds\n";
            flowCsv.open((outDir+"flow_stability.csv").c_str());
            flowCsv << "lb_step,t_s,max_velocity_lb,mean_velocity_lb,mean_density_w_lb,mean_density_o_lb,drive_factor,protocol_stage,min_total_density_lb,max_total_density_lb,inlet_outlet_density_rise_lb\n";
            protocolCsv.open((outDir+"restart_protocol.csv").c_str());
            protocolCsv << "lb_step,t_s,PV,stage,drive_factor,q_total_m3_s,core_dp_Pa,n_produced,n_weak_moving\n";
        }
    }
    std::ostringstream prodName; prodName << demOutDir << outputPrefix << "production_events_p" << global::mpi().getRank() << ".csv";
    ofstream prodCsv(prodName.str().c_str()); prodCsv << "lb_step,time_s,dem_step,id,type,x,y,z,radius\n";
    std::ostringstream satName; satName << demOutDir << outputPrefix << "particle_saturation_p" << global::mpi().getRank() << ".csv";
    if (calibrationMode==0) { satCsv.open(satName.str().c_str()); satCsv << "lb_step,t_s,id,type,x,local_Sw,immersion_time_s,n_sample_cells\n"; }

    bool outletOpen=false, completedPV=false, unstable=false, demNeedsSetup=false, calibrationConverged=false, calibrationAccepted=false;
    T finalCalibrationMd=(T)-1.;
    double producedWeak=0., cumulativeBrokenBonds=0., activeBonds=initialActiveBonds, weakNow=weakCount0;
    const plint equilSteps = calibrationMode==1 ? 0 : max<plint>(units.getLbSteps(equilibrationTime),1);
    T previousCalibrationK=-1., previousCalibrationQ=-1., previousCalibrationDp=-1.;
    std::vector<T> calibrationKWindow, calibrationQWindow, calibrationDpWindow; plint steadyCount=0;
    int lastStage=-2;
    bool densityWarned=false;      // ★★ G×ρ_总 过大只警告一次
    const T cs2_ = (T)1/(T)3;

    // ★★ 2026-10-26：各阶段耗时统计（墙钟时间，秒），每个日志步打印一次 [TIMING]，方便判断瓶颈
    double tSpheres=0., tShanChen=0., tCollide=0., tForces=0., tDem=0., tLogStart=MPI_Wtime();
    for (plint iT=0; iT<maxSteps; ++iT) {

        // (1) 颗粒 → 两套格子（固含率、颗粒速度）
        //     ★★ 2026-10-25：孔隙率场只在第 (8) 步"需要输出"时才刷新（旧版每步全域扫一遍）
        double tPhase = MPI_Wtime();
        if (calibrationMode==0) {
            setSpheresOnBothLattices(latticeW,latticeO,wrapper,units,&ibCouplingDomain,staticSolidFraction,staticParticleId);
        }
        tSpheres += MPI_Wtime()-tPhase;
        if (iT%vtkSteps==0) writeTwoPhaseVTK(latticeW,latticeO,units,iT,shanChenG,pScale,cs2_*rhoTotalBulk,outputPrefix);

        // (2) 出口回收：细粒球心越过 x=lx 即删除并记为产出
        if (calibrationMode==0 && !outletOpen && iT>=equilSteps && PV>=openOutletPV) { outletOpen=true; pcout << "[INFO] step=" << iT << " 出口回收区已启用\n"; }
        if (calibrationMode==0 && outletOpen && iT%particleCleanupInterval==0) {
            double producedLocal=0.;
            for (int ia=0; ia<wrapper.lmp->atom->nlocal; ++ia)
                if (wrapper.lmp->atom->type[ia]==2 && wrapper.lmp->atom->x[ia][0]>=lx) {
                    producedLocal += 1.;
                    prodCsv << iT << "," << iT*dt_phys << "," << wrapper.lmp->update->ntimestep << "," << wrapper.lmp->atom->tag[ia] << ",2,"
                            << wrapper.lmp->atom->x[ia][0] << "," << wrapper.lmp->atom->x[ia][1] << "," << wrapper.lmp->atom->x[ia][2] << "," << wrapper.lmp->atom->radius[ia] << "\n";
                }
            prodCsv.flush();
            double producedGlobal=0.; MPI_Allreduce(&producedLocal,&producedGlobal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            producedWeak += producedGlobal;
            if (producedGlobal>0.) {
                wrapper.execCommand("delete_atoms region prodzone compress no");
                demNeedsSetup=true;
                setSpheresOnBothLattices(latticeW,latticeO,wrapper,units,&ibCouplingDomain,staticSolidFraction,staticParticleId);
            }
        }

        // (3) 驱动因子（缓启 / 停注 / 分级重启）
        const T driveTime = max((T)0.,((T)(iT-equilSteps))*dt_phys);
        int stage=-1;
        const T factor = (calibrationMode==1) ? (T)1.0 : protocol.factor(driveTime, PV, stage);
        if (stage!=lastStage) { pcout << "[PROTOCOL] step=" << iT << " 进入阶段 " << stage << " 驱动因子=" << factor << "\n"; lastStage=stage; }

        // (4) 入口（恒流量或恒压力）与出口组成
        { vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&latticeO);
          applyProcessingFunctional(new OilWaterInlet3D(pressureControl?1:0, uInLb*factor, rhoTotalBulk+inletDeltaRho*factor, sInj, calibrationMode==1?true:waterInjected), inlet, a); }
        { vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&latticeO);
          applyProcessingFunctional(new OilWaterOutletComposition3D(rhoTotalBulk, rhoResidual), outlet, a); }

        // (5) 两相相互作用（公共速度 + 排斥力 + 润湿 + 毛细反作用力累加）
        tPhase = MPI_Wtime();
        if (calibrationMode==0) {
            capBuffers.prepare(wrapper, units);
            applyOilWaterShanChen(latticeW,latticeO,shanChenDomain,shanChenG,omegaW,omegaO,rhoResidual,wallK,containerWallK,
                                  mixtureTau==1, mixSatLow, mixSatHigh, traceVelocityCap, diffusionTau,   // ★★ 2026-10-25/26 新增参数
                                  &wrapper, capillaryOnParticles ? &capBuffers : 0,
                                  rhoWField, rhoOField, solidFlagField, scSolidThreshold);
        }

        tShanChen += MPI_Wtime()-tPhase;

        // (6) 碰撞 + 迁移（两套格子）
        tPhase = MPI_Wtime();
        latticeW.collideAndStream();
        latticeO.collideAndStream();
        tCollide += MPI_Wtime()-tPhase;
        if (iT < diagnosticScanSteps) diagnosticScan(latticeW, latticeO, iT);   // ★★ 排错：打印 [DIAG] 行

        // (7) 流量：分相统计入口截面的体积通量（动量 j_x = ρ u_x，油水格子密度都为 1，直接取平均）
        T qW_lb=0., qO_lb=0.;
        { auto rhoW = computeDensity(latticeW,fluxPlane); auto uxW = computeVelocityComponent(latticeW,fluxPlane,0);
          qW_lb = computeAverage(*multiply(*rhoW,*uxW));
          auto rhoO = computeDensity(latticeO,fluxPlane); auto uxO = computeVelocityComponent(latticeO,fluxPlane,0);
          qO_lb = computeAverage(*multiply(*rhoO,*uxO)); }
        if (!std::isfinite((double)qW_lb) || !std::isfinite((double)qO_lb)) { pcout << "[ERROR] step=" << iT << " 截面流量 NaN/Inf\n"; unstable=true; break; }
        const T qWater = qW_lb*units.getPhysVel(1)*channelAreaPhys, qOil = qO_lb*units.getPhysVel(1)*channelAreaPhys;
        const T volumeFlowRate = qWater+qOil;
        // ★★ 2026-10-26：恒流量模式下，注入量直接用入口给定的流量（V_INF × 驱动因子 × 截面积）累计 PV。
        //   旧版用 x=1 截面实测的注入相通量：开始驱替的头几百步里，初始界面松弛、两相溶解度调整引起的
        //   瞬态通量比注入速度大一两个数量级（本机测试剖面实测约 24 倍），PV 被高估，测试 70 步就"注满"了。
        //   恒压模式没有给定流量，仍用实测通量。实测通量照常写进 flow_stability.csv 供检查。
        const T injectedRate = pressureControl ? (waterInjected ? qWater : qOil)
                                               : (calibrationMode==1 ? volumeFlowRate : v_inf*factor*channelAreaPhys);
        if (iT>=equilSteps && injectedRate>0.) injectedVolume += injectedRate*dt_phys;
        PV = injectedVolume/initialPoreVolumePhys;

        // (8) 压降（★ P0-8：两相总压；标定模式退化为单相）
        //     ★★ 2026-10-25：只在"要写日志 / 要按 PV 采样 / 即将结束"的步才算。
        //     旧版每一步都新建两个全场标量场（N=40 时每个约 1.5 亿字节）并遍历全域，
        //     而结果只在这些步用到，白白多花约 20–40% 的时间。
        const bool pvSampleDue = (calibrationMode==0) && nextPV<pvSamples.size() && PV>=pvSamples[nextPV];
        const bool finishDue = (calibrationMode==0) && PV>=total_PV && !(protocol.stopRestart && stage<protocol.stages);
        const bool needFields = (iT%logSteps==0) || pvSampleDue || finishDue;
        T dpCore=-1., dpTotal=-1., screenDp=-1., mobility=-1., kAppMd=-1., densityRise=-1.;
        if (needFields) {
            if (calibrationMode==0) {   // 颗粒会动，先刷新孔隙率场（标定模式颗粒不动，用 t=0 的即可）
                vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&porosityField);
                applyProcessingFunctional(new ExtractPorosity3D(),latticeW.getBoundingBox(),a);
            }
            MultiScalarField3D<T> Sw(latticeW), pField(latticeW);
            { vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&latticeO); a.push_back(&Sw); a.push_back(&pField);
              applyProcessingFunctional(new OilWaterFields3D(shanChenG), latticeW.getBoundingBox(), a); }
            auto pWeighted = multiply(pField,porosityField);
            auto poreAvgP = [&](Box3D const& box)->T { T pore=computeSum(porosityField,box); return pore>0.?computeSum(*pWeighted,box)/pore:(T)0.; };
            dpCore = (poreAvgP(pressureLeft)-poreAvgP(pressureRight))*pScale;
            dpTotal = (poreAvgP(totalPressureIn)-poreAvgP(totalPressureOut))*pScale;
            screenDp = (poreAvgP(screenPressureUp)-poreAvgP(screenPressureDown))*pScale;
            mobility = (dpCore>1.e-15) ? (volumeFlowRate/channelAreaPhys)*pressureLengthPhys/dpCore : -1.;   // m²/(Pa·s)
            kAppMd = mobility>0. ? mobility*mu_f/(T)9.869233e-16 : (T)-1.;                                                       // 以水粘度折算的表观渗透率
            densityRise = dpTotal/pressurePerDensity;   // 实测入口−出口压差折成格子密度差（与预检的 Δρ 同一口径）
        }

        // (9) 标定模式：等单相渗透率收敛
        if (calibrationMode==1 && iT%logSteps==0 && iT>=equilSteps) {
            const T permeability = mobility*mu_f, permeabilityMd = kAppMd;
            T relativeChange=-1., windowRelativeRange=-1.;
            if (permeability>0. && previousCalibrationK>0. && volumeFlowRate>0. && previousCalibrationQ>0. && dpCore>0. && previousCalibrationDp>0.)
                relativeChange=max(fabs(permeability-previousCalibrationK)/previousCalibrationK,max(fabs(volumeFlowRate-previousCalibrationQ)/previousCalibrationQ,fabs(dpCore-previousCalibrationDp)/previousCalibrationDp));
            if (permeability>0. && volumeFlowRate>0. && dpCore>0.) {
                calibrationKWindow.push_back(permeability); calibrationQWindow.push_back(volumeFlowRate); calibrationDpWindow.push_back(dpCore);
                if ((plint)calibrationKWindow.size()>steadyWindowSize) { calibrationKWindow.erase(calibrationKWindow.begin()); calibrationQWindow.erase(calibrationQWindow.begin()); calibrationDpWindow.erase(calibrationDpWindow.begin()); }
            } else { calibrationKWindow.clear(); calibrationQWindow.clear(); calibrationDpWindow.clear(); }
            if ((plint)calibrationKWindow.size()==steadyWindowSize) {
                T minK=calibrationKWindow[0],maxK=minK,sumK=0.,minQ=calibrationQWindow[0],maxQ=minQ,sumQ=0.,minDp=calibrationDpWindow[0],maxDp=minDp,sumDp=0.;
                for (plint iw=0;iw<steadyWindowSize;++iw) {
                    minK=min(minK,calibrationKWindow[iw]);maxK=max(maxK,calibrationKWindow[iw]);sumK+=calibrationKWindow[iw];
                    minQ=min(minQ,calibrationQWindow[iw]);maxQ=max(maxQ,calibrationQWindow[iw]);sumQ+=calibrationQWindow[iw];
                    minDp=min(minDp,calibrationDpWindow[iw]);maxDp=max(maxDp,calibrationDpWindow[iw]);sumDp+=calibrationDpWindow[iw];
                }
                windowRelativeRange=max((maxK-minK)/(sumK/steadyWindowSize),max((maxQ-minQ)/(sumQ/steadyWindowSize),(maxDp-minDp)/(sumDp/steadyWindowSize)));
            }
            if (windowRelativeRange>=0. && windowRelativeRange<steadyRelativeTolerance) ++steadyCount; else steadyCount=0;
            previousCalibrationK=permeability; previousCalibrationQ=volumeFlowRate; previousCalibrationDp=dpCore;
            if (global::mpi().getRank()==0) {
                calibrationCsv << iT << "," << (T)iT*dt_phys << "," << phi0 << "," << volumeFlowRate << "," << dpCore << "," << dpTotal << "," << permeability << "," << permeabilityMd << "," << relativeChange << "," << windowRelativeRange << "," << steadyCount << "\n";
                calibrationCsv.flush();
            }
            pcout << "[PERMEABILITY_CALIBRATION] step=" << iT << " k=" << permeabilityMd << "mD coreDp=" << dpCore << "Pa steady=" << steadyCount << "/" << requiredSteadyCount << std::endl;
            if (steadyCount>=requiredSteadyCount) {
                calibrationConverged=true; finalCalibrationMd=permeabilityMd; calibrationAccepted=permeabilityMd>=initialKMinMd && permeabilityMd<=initialKMaxMd; completedPV=true;
                pcout << "[PERMEABILITY_RESULT] converged=1 permeability_md=" << permeabilityMd << " accepted=" << (calibrationAccepted?1:0) << std::endl;
                break;
            }
        }

        // (10) 生产模式：两相合力 → 力审计 → DEM 推进 → 守恒检查 → 颗粒饱和度采样
        if (calibrationMode==0) {
            T maxCapForce=0.;
            tPhase = MPI_Wtime();
            getTwoPhaseForcesFromLattice(latticeW,latticeO,wrapper,units,capillaryOnParticles?&capBuffers:0,maxCapForce);
            tForces += MPI_Wtime()-tPhase;
            T stepMaxCapForce=0.; MPI_Allreduce(&maxCapForce,&stepMaxCapForce,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

            LAMMPS_NS::FixLbCouplingOnetoone* couplingFix = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>(wrapper.lmp->modify->find_fix_style("couple/lb/onetoone",0));
            double** fL = couplingFix->get_force_ptr(); double** tL = couplingFix->get_torque_ptr();
            // ★★ 2026-10-26：受力审计只对"会动的细粒"（type 2）判阈值。
            //   旧版把固定骨架（type 1，直径约 0.3 mm，体积是细粒的 200 多倍）和筛网（type 3）也算进去，
            //   骨架颗粒在正常压力梯度/启动压力波下受力就可能超过为细粒设定的 IBM_FORCE_ABORT_N，
            //   而它们根本不动、受力也不进入 DEM（fix lbcoupling 只作用于 weak 组），却会让计算被误判中止。
            //   现在：细粒超阈值 → 中止；任何颗粒出现 NaN/Inf → 中止；骨架/筛网的最大受力只记录和打印。
            //   失败时额外打印受力最大的细粒编号、位置（格子坐标）和半径，方便定位。
            plint nPart = wrapper.lmp->atom->nlocal+wrapper.lmp->atom->nghost;
            T localMaxF=0.,localMaxT=0.,localMaxFixedF=0.; plint localNonFinite=0;
            double localTop[6]={0.,-1.,0.,0.,0.,0.};   // |F|, id, x, y, z (格子单位), 半径(格)
            double localTopFixed[6]={0.,-1.,0.,0.,0.,0.};   // 同上，固定骨架/筛网里受力最大的
            for (plint ip=0; ip<nPart; ++ip) {
                const bool mobile = wrapper.lmp->atom->type[ip]==2;
                T fMag=0.;
                for (int j=0;j<3;++j) {
                    double fv=fL[ip][j],tv=tL[ip][j];
                    if (!std::isfinite(fv)||!std::isfinite(tv)) ++localNonFinite;
                    fMag += fv*fv;
                    if (mobile) localMaxT=max(localMaxT,(T)fabs(tv));
                }
                fMag = std::sqrt(fMag);
                if (!mobile) {
                    localMaxFixedF=max(localMaxFixedF,fMag);
                    if (ip<wrapper.lmp->atom->nlocal && fMag>localTopFixed[0]) {
                        localTopFixed[0]=fMag; localTopFixed[1]=wrapper.lmp->atom->tag[ip]+0.1*wrapper.lmp->atom->type[ip];
                        for (int d=0;d<3;++d) localTopFixed[2+d]=wrapper.lmp->atom->x[ip][d]/physDx;
                        localTopFixed[5]=wrapper.lmp->atom->radius[ip]/physDx;
                    }
                    continue;
                }
                localMaxF=max(localMaxF,fMag);
                if (ip<wrapper.lmp->atom->nlocal && fMag>localTop[0]) {
                    localTop[0]=fMag; localTop[1]=wrapper.lmp->atom->tag[ip];
                    for (int d=0;d<3;++d) localTop[2+d]=wrapper.lmp->atom->x[ip][d]/physDx;
                    localTop[5]=wrapper.lmp->atom->radius[ip]/physDx;
                }
            }
            T stepMaxForce=0.,stepMaxTorque=0.,stepMaxFixedForce=0.; int stepNonFinite=0; int localNonFiniteInt=(int)localNonFinite;
            MPI_Allreduce(&localMaxF,&stepMaxForce,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            MPI_Allreduce(&localMaxT,&stepMaxTorque,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            MPI_Allreduce(&localMaxFixedF,&stepMaxFixedForce,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            MPI_Allreduce(&localNonFiniteInt,&stepNonFinite,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            bool overForce = stepMaxForce>ibmForceAbort, overTorque = stepMaxTorque>ibmTorqueAbort;
            if (stepNonFinite>0 || overForce || overTorque) {
                pcout << "[ERROR] step=" << iT << " 受力审计失败: 细粒 maxF=" << stepMaxForce << " N  maxT=" << stepMaxTorque << " N·m  nonFinite=" << stepNonFinite
                      << "（阈值 IBM_FORCE_ABORT_N=" << ibmForceAbort << " / IBM_TORQUE_ABORT_NM=" << ibmTorqueAbort << "）；固定骨架/筛网最大受力=" << stepMaxFixedForce << " N\n";
                {   // 各进程把自己受力最大的细粒交给 0 号进程打印
                    int nProc=global::mpi().getSize();
                    std::vector<double> all(6*nProc,0.), allFixed(6*nProc,0.);
                    MPI_Gather(localTop,6,MPI_DOUBLE,&all[0],6,MPI_DOUBLE,0,MPI_COMM_WORLD);
                    MPI_Gather(localTopFixed,6,MPI_DOUBLE,&allFixed[0],6,MPI_DOUBLE,0,MPI_COMM_WORLD);
                    if (global::mpi().getRank()==0) {
                        int best=0; for (int r=1;r<nProc;++r) if (all[6*r]>all[6*best]) best=r;
                        if (all[6*best+1]>0)
                            std::cout << "        受力最大的细粒: id=" << (long)all[6*best+1] << " |F|=" << all[6*best] << " N  位置(格)=("
                                      << all[6*best+2] << "," << all[6*best+3] << "," << all[6*best+4] << ")  半径=" << all[6*best+5] << " 格" << std::endl;
                        best=0; for (int r=1;r<nProc;++r) if (allFixed[6*r]>allFixed[6*best]) best=r;
                        if (allFixed[6*best+1]>0) {
                            long const idType = (long)std::floor(allFixed[6*best+1]);
                            int const ty = (int)std::floor((allFixed[6*best+1]-idType)*10.+0.5);
                            std::cout << "        受力最大的固定颗粒: id=" << idType << " type=" << ty << " |F|=" << allFixed[6*best] << " N  位置(格)=("
                                      << allFixed[6*best+2] << "," << allFixed[6*best+3] << "," << allFixed[6*best+4] << ")  半径=" << allFixed[6*best+5] << " 格" << std::endl;
                        }
                    }
                }
                unstable=true;
                // ★★ 2026-10-25：旧版这里只做标记、不停算，DEM 被冻结而 LBM 继续跑到注满 PV，最后还以 0 退出。
                //    现在写一行审计记录、存最后一帧 VTK 后立即停算（退出码 2）。
                if (global::mpi().getRank()==0) {
                    forceCsv << iT << "," << (T)iT*dt_phys << "," << stepMaxForce << "," << stepMaxTorque << "," << stepMaxCapForce << ",-1,-1,-1,"
                             << (overForce?1:0) << "," << (overTorque?1:0) << "," << stepNonFinite << "," << stepMaxFixedForce << "\n";
                    forceCsv.flush();
                }
                writeTwoPhaseVTK(latticeW,latticeO,units,iT,shanChenG,pScale,cs2_*rhoTotalBulk,outputPrefix);
                break;
            } else {
                tPhase = MPI_Wtime();
                if (demNeedsSetup) { wrapper.run(demSubsteps); demNeedsSetup=false; } else wrapper.runWithoutSetup(demSubsteps);
                tDem += MPI_Wtime()-tPhase;
                bondCompute->compute_vector();
                double brokenThisStep=global::mpi().getRank()==0?bondCompute->vector[1]:0.;
                activeBonds=global::mpi().getRank()==0?bondCompute->vector[2]:0.;
                MPI_Bcast(&brokenThisStep,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&activeBonds,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                cumulativeBrokenBonds += brokenThisStep;
                if (fabs(activeBonds+cumulativeBrokenBonds-initialActiveBonds)>0.1) { pcout << "[ERROR] step=" << iT << " 胶结键计数不守恒\n"; unstable=true; }
                T localMaxParticleSpeed=0.,localMaxParticleSurfaceSpeed=0.; double movingLocal=0.;
                for (int ia=0;ia<wrapper.lmp->atom->nlocal;++ia) {
                    if (wrapper.lmp->atom->type[ia]!=2) continue;
                    const T vx=wrapper.lmp->atom->v[ia][0],vy=wrapper.lmp->atom->v[ia][1],vz=wrapper.lmp->atom->v[ia][2];
                    const T wx=wrapper.lmp->atom->omega[ia][0],wy=wrapper.lmp->atom->omega[ia][1],wz=wrapper.lmp->atom->omega[ia][2];
                    const T sp=sqrt(vx*vx+vy*vy+vz*vz);
                    localMaxParticleSpeed=max(localMaxParticleSpeed,sp);
                    localMaxParticleSurfaceSpeed=max(localMaxParticleSurfaceSpeed,wrapper.lmp->atom->radius[ia]*sqrt(wx*wx+wy*wy+wz*wz));
                    if (sp>screenBlockedSpeed) movingLocal+=1.;
                }
                T stepMaxParticleSpeed=0.,stepMaxParticleSurfaceSpeed=0.; double movingNow=0.;
                MPI_Allreduce(&localMaxParticleSpeed,&stepMaxParticleSpeed,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                MPI_Allreduce(&localMaxParticleSurfaceSpeed,&stepMaxParticleSurfaceSpeed,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                MPI_Allreduce(&movingLocal,&movingNow,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                T velocityLimitCount=integratorFix->compute_scalar();
                if (stepMaxParticleSpeed>demMaxParticleSpeed*1.000001 || stepMaxParticleSurfaceSpeed>demMaxParticleSpeed*1.000001) { pcout << "[ERROR] step=" << iT << " DEM速度限制失效: v=" << stepMaxParticleSpeed << "\n"; unstable=true; }

                double weakAfterLocal=0.,skeletonAfterLocal=0.,screenAfterLocal=0.;
                for (int ia=0;ia<wrapper.lmp->atom->nlocal;++ia) {
                    if (wrapper.lmp->atom->type[ia]==2) weakAfterLocal+=1.;
                    if (wrapper.lmp->atom->type[ia]==1) skeletonAfterLocal+=1.;
                    if (wrapper.lmp->atom->type[ia]==3) screenAfterLocal+=1.;
                }
                double skeletonNow=0.,screenNow=0.;
                MPI_Allreduce(&weakAfterLocal,&weakNow,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&skeletonAfterLocal,&skeletonNow,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&screenAfterLocal,&screenNow,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                const double abnormalAfter=weakCount0-producedWeak-weakNow;
                if (fabs(abnormalAfter)>0.1 || fabs(skeletonNow-skeletonCount0)>0.1 || fabs(screenNow-screenCount0)>0.1) { pcout << "[ERROR] step=" << iT << " DEM颗粒守恒失败\n"; unstable=true; }

                // ★ 颗粒周围含水饱和度采样（每个 dump 间隔一次），写到每进程 CSV
                if (iT%particleDumpSteps==0) {
                    satSumW.assign(nPart,(T)0); satSumO.assign(nPart,(T)0);
                    plint const iBlock = latticeW.getLocalInfo().getBlocks()[0];
                    std::map<plint,Box3D> blockmap = latticeW.getSparseBlockStructure().getBulks();
                    vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&latticeO);
                    applyProcessingFunctional(new ParticleSaturationSampler3D(&satSumW,&satSumO,&wrapper), blockmap[iBlock], a);
                    for (int ia=0; ia<wrapper.lmp->atom->nlocal; ++ia) {
                        if (wrapper.lmp->atom->type[ia]!=2) continue;
                        T const tot = satSumW[ia]+satSumO[ia];
                        T const swLocal = tot>0. ? satSumW[ia]/tot : -1.;
                        int const tag = wrapper.lmp->atom->tag[ia];
                        if (swLocal>0.5) immersionTime[tag] += particleDumpSteps*dt_phys;
                        satCsv << iT << "," << iT*dt_phys << "," << tag << ",2," << wrapper.lmp->atom->x[ia][0] << "," << swLocal << "," << immersionTime[tag] << "," << tot << "\n";
                    }
                    satCsv.flush();
                }

                if (global::mpi().getRank()==0 && (iT%logSteps==0 || unstable || PV>=total_PV)) {
                    forceCsv << iT << "," << (T)iT*dt_phys << "," << stepMaxForce << "," << stepMaxTorque << "," << stepMaxCapForce << "," << stepMaxParticleSpeed << "," << stepMaxParticleSurfaceSpeed << "," << velocityLimitCount << "," << (overForce?1:0) << "," << (overTorque?1:0) << "," << stepNonFinite << "," << stepMaxFixedForce << "\n"; forceCsv.flush();
                    balanceCsv << iT << "," << (T)iT*dt_phys << "," << skeletonCount0 << "," << skeletonNow << "," << screenCount0 << "," << screenNow << "," << weakCount0 << "," << weakNow << "," << producedWeak << "," << abnormalAfter << "\n"; balanceCsv.flush();
                    bondCsv << iT << "," << (T)iT*dt_phys << "," << initialActiveBonds << "," << activeBonds << "," << brokenThisStep << "," << cumulativeBrokenBonds << "\n"; bondCsv.flush();
                    protocolCsv << iT << "," << (T)iT*dt_phys << "," << PV << "," << stage << "," << factor << "," << volumeFlowRate << "," << dpCore << "," << producedWeak << "," << movingNow << "\n"; protocolCsv.flush();
                }
                if (unstable) break;
            }

            // (11) 按 PV 采样孔渗
            while(nextPV<pvSamples.size() && PV>=pvSamples[nextPV]) {
                const T phi_now=computeSum(porosityField,matrixCore)/matrixGeometricNodes;
                pcout << "[PV] targetPV=" << pvSamples[nextPV] << " actualPV=" << PV << " mobility=" << mobility << " k_app(μ_w)=" << kAppMd << "mD coreDp=" << dpCore << " totalDp=" << dpTotal << " nProd=" << producedWeak << std::endl;
                if (global::mpi().getRank()==0) {
                    permCsv << PV << "," << pvSamples[nextPV] << "," << (T)iT*dt_phys << "," << phi_now << "," << mobility << "," << kAppMd << "," << dpCore << "," << dpTotal << "," << qWater << "," << qOil << "," << producedWeak << "," << weakNow << "\n"; permCsv.flush();
                }
                nextPV++;
            }
        }

        if (calibrationMode==0 && PV>=total_PV && !(protocol.stopRestart && stage<protocol.stages)) {
            pcout << "[PV] 达到总注入 " << total_PV << " PV (step " << iT << ")\n";
            writeTwoPhaseVTK(latticeW,latticeO,units,iT,shanChenG,pScale,cs2_*rhoTotalBulk,outputPrefix);
            completedPV=true; break;
        }

        // (12) 日志：筛网压降、流场稳定性、密度监测
        if (iT%logSteps==0) {
            double nearScreenLocal=0.,nearScreenSlowLocal=0.,downstreamScreenLocal=0.;
            for (int ia=0;ia<wrapper.lmp->atom->nlocal;++ia) {
                if (wrapper.lmp->atom->type[ia]!=2) continue;
                const T x=wrapper.lmp->atom->x[ia][0],radius=wrapper.lmp->atom->radius[ia];
                const T vx=wrapper.lmp->atom->v[ia][0],vy=wrapper.lmp->atom->v[ia][1],vz=wrapper.lmp->atom->v[ia][2];
                const T speed=sqrt(vx*vx+vy*vy+vz*vz);
                if (x<screenX && x+radius>=screenX-screenCaptureBand) { nearScreenLocal+=1.; if (speed<screenBlockedSpeed) nearScreenSlowLocal+=1.; }
                if (x-radius>screenX) downstreamScreenLocal+=1.;
            }
            double nearScreen=0.,nearScreenSlow=0.,downstreamScreen=0.;
            MPI_Allreduce(&nearScreenLocal,&nearScreen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(&nearScreenSlowLocal,&nearScreenSlow,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(&downstreamScreenLocal,&downstreamScreen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            if (global::mpi().getRank()==0) { screenCsv << iT << "," << (T)iT*dt_phys << "," << PV << "," << screenDp << "," << dpTotal << "," << nearScreen << "," << nearScreenSlow << "," << downstreamScreen << "," << producedWeak << "\n"; screenCsv.flush(); }
            if (calibrationMode==0) {
                // ★★ 2026-10-25：用两相合起来的速度判失稳（见 TwoPhaseStats3D 的说明），并监测总密度
                TwoPhaseStats3D stats;
                applyProcessingFunctional(stats, shanChenDomain, latticeW, latticeO);
                T const uMaxLb = stats.maxSpeed(), uMeanLb = stats.meanSpeed();
                T const rhoMeanW = stats.meanRhoW(), rhoMeanO = stats.meanRhoO();
                T const rhoTotMin = stats.minRho(), rhoTotMax = stats.maxRho();
                if (global::mpi().getRank()==0) {
                    flowCsv << iT << "," << (T)iT*dt_phys << "," << uMaxLb << "," << uMeanLb << "," << rhoMeanW << "," << rhoMeanO << "," << factor << "," << stage
                            << "," << rhoTotMin << "," << rhoTotMax << "," << densityRise << "\n";
                    flowCsv.flush();
                }
                pcout << "[DENSITY] step=" << iT << " 入口−出口实测 Δρ=" << densityRise << "（预检估计 " << predictedDensityDrop
                      << "）  孔隙内 ρ_总 ∈ [" << rhoTotMin << ", " << rhoTotMax << "]  G×ρ_总,max=" << shanChenG*rhoTotMax << "\n";
                if (stats.nanCount()>0 || !(uMaxLb<0.35)) {
                    pcout << "[ERROR] step=" << iT << " 流场失稳: 混合物 max|u|=" << uMaxLb << "（格子单位，应 <0.35）, NaN 格点数=" << stats.nanCount() << "\n";
                    writeTwoPhaseVTK(latticeW,latticeO,units,iT,shanChenG,pScale,cs2_*rhoTotalBulk,outputPrefix);
                    unstable=true; break;
                }
                if (densityRise > maxRuntimeDensityRise) {
                    pcout << "[ERROR] step=" << iT << " 实测入口−出口格子密度差 Δρ=" << densityRise << " 超过 MAX_RUNTIME_LB_DENSITY_RISE="
                          << maxRuntimeDensityRise << "：压差太大，Shan-Chen 的界面张力已沿岩心明显变化、接近发散。\n"
                          << "        通常是因为实际渗透率低于 K_MIN_MD 或细粒堵塞筛网。请降低 V_INF（恒流）或 INLET_PRESSURE_PA（恒压）。\n";
                    writeTwoPhaseVTK(latticeW,latticeO,units,iT,shanChenG,pScale,cs2_*rhoTotalBulk,outputPrefix);
                    unstable=true; break;
                }
                if (shanChenG*rhoTotMax > (T)3.3 && !densityWarned) {
                    pcout << "[WARN] G×ρ_总 最大值=" << shanChenG*rhoTotMax << " > 3.3（G=3.5 在验证中发散），请留意后续是否失稳。\n";
                    densityWarned = true;
                }
            }
            // ★★ 2026-10-26：改用墙钟时间（旧版 clock() 是 CPU 时间，多进程忙等时会偏大）
            const double now = MPI_Wtime(); const double elapsed = now - tLogStart;
            const double nSteps = (iT==0) ? 1. : (double)logSteps;
            pcout << "[step=" << iT << "] t=" << iT*dt_phys << "s PV=" << PV << " " << elapsed << "s "
                  << ((double)(latticeW.getNx()*latticeW.getNy()*latticeW.getNz())*nSteps)/elapsed/1e6 << " MLU/s\n";
            pcout << "[TIMING] 每步平均耗时(s)：颗粒映射=" << tSpheres/nSteps << " 两相内核=" << tShanChen/nSteps
                  << " 碰撞迁移=" << tCollide/nSteps << " 颗粒受力=" << tForces/nSteps << " DEM=" << tDem/nSteps
                  << " 合计=" << elapsed/nSteps << "\n";
            tSpheres=tShanChen=tCollide=tForces=tDem=0.; tLogStart=MPI_Wtime();
            {   // ★★ 2026-10-27：痕量组分密度地板修复次数（见 oilwaterShanChen3D.hh 修订八）
                long nLoc = traceDensityRepairCount(), nGlob = 0; traceDensityRepairCount() = 0;
                MPI_Allreduce(&nLoc,&nGlob,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
                if (nGlob > 0) {
                    const double perStep = (double)nGlob/nSteps;
                    const double ppm = perStep/(double)(latticeW.getNx()*latticeW.getNy()*latticeW.getNz())*1e6;
                    pcout << "[REPAIR] 本日志段痕量组分密度修复 " << nGlob << " 次，平均每步 " << perStep << " 个格点（占全场 "
                          << ppm << " ppm）。几到几十 ppm 属正常（每次只补回 <0.005 的质量）；若持续 >1000 ppm 且逐段增长，"
                          << "说明界面附近数值不稳，请把日志发给我\n";
                }
            }
        }
    }

    // ── 收尾 ─────────────────────────────────────────────────────────────
    if (global::mpi().getRank()==0) {
        ofstream* files[] = {&permCsv,&calibrationCsv,&forceCsv,&balanceCsv,&bondCsv,&screenCsv,&flowCsv,&protocolCsv};
        for (int i=0;i<8;++i) if (files[i]->is_open()) files[i]->close();
    }
    if (prodCsv.is_open()) prodCsv.close();
    if (satCsv.is_open()) satCsv.close();
    { std::stringstream finalDump; finalDump << "write_dump all custom " << demOutDir << outputPrefix << "particles_final.liggghts id type x y z vx vy vz radius fx fy fz modify sort id"; wrapper.execCommand(finalDump); }

    int conversionFailed=0;
    if (global::mpi().getRank()==0) {
        std::string cmd = "for f in \"" + demOutDir + "\"" + outputPrefix + "particles_*.liggghts; do python3 convert_liggghts_to_csv.py \"$f\"; done";
        if (system(cmd.c_str())!=0) conversionFailed=1;
    }
    MPI_Bcast(&conversionFailed,1,MPI_INT,0,MPI_COMM_WORLD);
    if (unstable && completedPV) { pcout << "[ERROR] 运行中出现过失稳/审计失败，结果不可信\n"; return 2; }
    if (!completedPV) {
        if (!unstable && calibrationMode==1) pcout << "[ERROR] 达到 maxT 但渗透率未收敛\n";
        else if (!unstable) pcout << "[ERROR] 达到 maxT 但仅注入 " << PV << "/" << total_PV << " PV\n";
        return unstable?2:3;
    }
    if (calibrationMode==1) {
        if(conversionFailed) return 5;
        if(!calibrationConverged) return 3;
        if(!calibrationAccepted) { pcout << "[ERROR] 稳态渗透率=" << finalCalibrationMd << "mD 不在区间 " << initialKMinMd << ".." << initialKMaxMd << "\n"; return 4; }
        return 0;
    }
    if (conversionFailed) return 5;
    int postFailed=0;
    if (global::mpi().getRank()==0) { std::string cmd2 = "python3 analyze_production.py \"" + outDir + "\""; if (system(cmd2.c_str())!=0) postFailed=1; }
    MPI_Bcast(&postFailed,1,MPI_INT,0,MPI_COMM_WORLD);
    if (postFailed) return 6;
    delete bcW; delete bcO;
    return 0;
}

/*
 * oilwaterShanChen3D.h
 * ====================================================================================
 * 油–水两相 Shan-Chen 多组分相互作用处理器（修订版）
 *
 * 【给初学者的说明】
 *   程序里油和水各有一套格子（latticeW、latticeO），每套格子每个格点存 19 个方向的
 *   分布函数 f_i。Shan-Chen 模型用一个"跨组分排斥力"让油水自动分开：
 *       F_w(x) = −G · ρ_w(x) · Σ_i t_i ρ_o(x+c_i) c_i      （水被油"推开"）
 *       F_o(x) = −G · ρ_o(x) · Σ_i t_i ρ_w(x+c_i) c_i      （油被水"推开"）
 *   G 越大排斥越强、界面张力越大；G 太大会发散。
 *
 * 【与原版（审计前）的三处本质区别】——这三条正是审计报告 P0-1、P0-4、P0-5
 *   1. 公共速度 u'：油、水碰撞时都向同一个速度
 *          u' = (ω_w j_w + ω_o j_o) / (ω_w ρ_w + ω_o ρ_o)
 *      靠拢（Shan & Doolen 1995）。原版各用各的速度，两相之间没有粘性耦合，
 *      按原参数第 1 步就失稳（见审计附录 A 算例 1）。
 *   2. 精确的力注入：原版把 τ·F 写进分布函数，每步实际注入 τF 而不是 F，
 *      界面张力会随 τ 变、油水不对称。本版把碰撞前的动量调到
 *          j* = (1−ω) j + ω ρ u' + F
 *      ——这正好是 Shan-Doolen BGK 碰撞后应有的动量；后面的正则化碰撞保持动量，
 *      所以每步净注入恰好是 F。
 *   3. 润湿性（接触角）：固体格点（颗粒/壁面内部）不再用它自己漂移的分布函数，
 *      而用一个人为规定的"虚拟密度"参与求和。原版固体格点里存的是初始化时"碰巧"填进去的
 *      油或水，所以岩心默认是强油湿且无法控制（审计附录 A 算例 3）。
 *      ★★ 2026-10-26 起虚拟密度按"几何接触角"给出（修订六，见下），参数直接是接触角
 *      WALL_CONTACT_ANGLE_DEG_TYPE1/2/3（从水相一侧量，<90° 水湿，>90° 油湿）。
 *
 * 【2026-10-25 修订一：局部混合松弛时间（解决粘度比 M≳9 时发散）】
 *   原做法：水格子用 τ_w、油格子用 τ_o，各自固定。验证发现有壁面时 τ_o≥4.5（τ_w=0.8 时 M≳10）
 *   在 20 步内发散——油格子松弛太慢，残余组分摆脱公共速度后失控。
 *   现做法（mixTau=true，参数 MIXTURE_TAU=1）：每个格点按局部含水比例 s=ρ_w/(ρ_w+ρ_o)
 *   求一个"混合粘度"，两套格子在该格点用同一个 ω：
 *       φ = clamp((s − s_lo)/(s_hi − s_lo), 0, 1)       （s_lo/s_hi = MIX_SAT_LOW/HIGH）
 *       1/ν_mix = φ/ν_w + (1−φ)/ν_o                      （调和平均，界面处平滑过渡）
 *       ω_mix = 1/(0.5 + 3 ν_mix)
 *   本体水区 ν=ν_w、本体油区 ν=ν_o，粘度比不变；numpy 验证（validation/v5_high_viscosity_ratio.py）
 *   在 τ_o=19（M≈50）下仍稳定。实现方式：每步对每个流体格点调用 setOmega(ω_mix)，
 *   PSM 的 B 加权也随之使用同一个 ω，两相看到的颗粒边界从此一致。
 *
 * 【2026-10-25 修订二：痕量组分速度上限只限制"偏离公共速度"的部分】
 *   旧做法把痕量组分的整个目标速度 j_star/ρ 按比例缩到 0.3 以内。界面两侧的痕量组分
 *   受到很强的法向排斥（F/ρ ≈ G·c_s²·|∇ρ_另一相| ≈ 0.3~0.4），一旦被缩放，
 *   它沿界面方向（切向）的速度也被一起缩小 → 每步都在界面上"吃掉"切向动量，
 *   界面变成一堵"无滑移墙"。numpy 分层流验证：同粘度（M=1）两层流误差 96%，
 *   界面处速度几乎为 0。
 *   现做法：只限制 (j_star/ρ − u') 这个"偏离公共速度"的部分，公共速度 u' 原样保留：
 *       v = j_star/ρ − u' ；若 |v| > uCap，则 v ← v·uCap/|v| ；j_star = ρ (u' + v)
 *   法向的稳定作用不变（壁面润湿强时仍需要它），切向不再损失动量。改后 M=1 分层流
 *   L2 误差 4.5%、流量误差 1.7%。上限值由参数 TRACE_VELOCITY_CAP 给出（默认 0.3）。
 *   注意：这一改法只在 MIXTURE_TAU=1 时使用。MIXTURE_TAU=0（各用固定 τ）时保留旧的整体缩放——
 *   固定 τ 的做法本身在壁面附近不够稳定，旧缩放的额外阻尼正是它能算下去的原因
 *   （固定 τ + 只限偏离量，在 s_wall≥0.5 的接触角算例里发散）。
 *
 * 【2026-10-26 修订五：组分间动量交换用单独的松弛时间 τ_D（参数 DIFFUSION_TAU，默认 1.0）】
 *   旧做法里"公共速度 u'"和"j* = (1−ω) j + ω ρ u' + F"用的 ω 就是粘性 ω（混合 τ）。油区 τ≈5.45，
 *   ω≈0.18，痕量组分偏离公共速度的稳态值 F/(ωρ) 是 τ=1 时的 5.5 倍，经常顶到 TRACE_VELOCITY_CAP 被削掉，
 *   削掉的动量不守恒 → 界面处出现净"泵送力"。numpy 一维平界面验证（τ 0.7/5.45，2 万步）：
 *   油、水两侧压力差 0.019（格子单位）、整体漂移速度 0.012（注入速度的 1000 多倍）。
 *   参照 Porter 等（2012, PRE 86:036701）"组分扩散与粘度分开松弛"的思路：动量交换（决定组分间扩散）
 *   统一用 ω_D = 1/τ_D，粘性应力仍按混合 τ 在随后的正则化碰撞里松弛。这样两相动量严格守恒，
 *   一维平界面两侧压力相等、残余速度 2.5e-6；粘度比不受影响（分层 Poiseuille 验证见 validation/）。
 *   DIFFUSION_TAU=0 可切回旧做法做对比。
 *
 * 【2026-10-26 修订六：几何接触角润湿模型】
 *   旧的两种虚拟密度写法都会在整个固体表面对"本体"流体施力：
 *     · ρ̃_w = s·ρ_wall、ρ̃_o = (1−s)·ρ_wall：壁面同时排斥两相，窄喉道被挤爆（本机第 7 步发散）；
 *     · 修订三"只排斥非润湿相"：中性壁面没问题，但油湿/水湿时非润湿相整体被推离壁面、在喉道里被压缩，
 *       numpy 二维多孔介质算例 s=0.3、0.7 都在 500 步内发散。
 *   现做法（Ding & Spelt 2007；Li 等 2014 用于伪势模型）：相场 φ = (ρ_w−ρ_o)/(ρ_w+ρ_o)，
 *   壁面法向 n 由固体邻格方向求出，固体邻格 x+c_i 的 φ 按接触角条件外推：
 *       φ_s = φ(x) + c_i·∇_tφ − K·|∇_tφ|·(c_i·n)，K = WETTING_GAIN·cotθ
 *   总密度取本格的值（镜像）。本体流体里 ∇φ≈0，壁面不施任何附加力；只有接触线附近起作用。
 *   numpy 验证（validation/v3_contact_angle.py，τ 0.7/5.45）：输入角与实测角单调对应，
 *   二维多孔介质算例在 θ=30°～150° 全部稳定（旧写法 s=0.3/0.7 都发散）。
 *
 * 【2026-10-26 修订七：Shan-Chen 眼里的"固体"以固含率 0.5 为界（参数 SHAN_CHEN_SOLID_THRESHOLD）】
 *   旧做法只把"完全被颗粒覆盖"（固含率 > 0.999）的格点当固体，固含率 0.9、0.99 的格点仍当流体，
 *   照样施加跨组分力和润湿。这些格点里的"流体"几乎全在颗粒里，PSM 每步都把它们的动量拉回颗粒速度；
 *   在细粒与骨架相切处的楔形缝里，一个格点周围同时有两个颗粒的高固含率格点，Shan-Chen 力和 PSM 反弹
 *   互相"打架"（本机测试剖面第 2159 步正是在细粒 321 与骨架 100 的切点处发散）。
 *   现做法：固含率 > 0.5 的格点在 Shan-Chen 里一律当固体（参与几何接触角的"虚拟密度"），不再对它施力；
 *   颗粒表面相当于取在固含率 0.5 的等值面上（PSM 的常用约定）。这些格点仍按局部含水比例设置混合粘度，
 *   保证 PSM 的 B 加权和流体碰撞用的 τ 与两相一致。
 *
 * 【2026-10-27 修订八：痕量组分的松弛与"密度地板"】
 *   现象：本机测试剖面第 1338–1343 步，骨架颗粒 63 旁（固含率 0.18 的格点）痕量水的密度先在 −0.0005 左右
 *   停留上百步，随后 0.002 → 0.49 → 15.5 → −5.6e5 爆掉。
 *   原因：① 油区里的痕量水也按油的粘性 τ≈5.45 松弛，它的非平衡部分（应力）衰减极慢，PSM 反弹和几何润湿
 *   每步给它一点扰动，就会累积成负的分布函数和量；② 密度一旦为负，Shan-Chen 力和公共速度都失去意义，
 *   这个格点之后只会越来越坏。
 *   修正：(a) 某组分的质量分数低于 MIX_SAT_LOW 时，它的松弛频率从混合 ω 逐渐过渡到 ω_D = 1/τ_D（τ_D=1）：
 *            ω_组分 = w·ω_mix + (1−w)·ω_D，  w = min(1, 该组分质量分数 / MIX_SAT_LOW)
 *         本体组分（w=1）完全不变，所以粘度比、分层流验证结果不受影响；痕量组分本来就只承担"跟着走"的角色，
 *         让它快些回到平衡态即可（与修订五的 τ_D 同一思路）。
 *       (b) 痕量组分密度低于 0.5·RHO_RESIDUAL 时，把它重置为"密度 RHO_RESIDUAL、速度=另一组分速度"的平衡态
 *         （流体格点和颗粒覆盖的格点都做）。每次重置会给这个组分补回 < 0.005 的质量，日志 [REPAIR] 行
 *         打印累计次数；正常计算里这个数应远小于格点数×步数（详见 修改说明 补充三）。
 *
 * 【新增：毛细力对颗粒的反作用】——审计 P0-4
 *   流体格点 x 受到固体邻格（属于颗粒 p）的虚拟密度产生的排斥力 dF，按牛顿第三定律，
 *   颗粒 p 受到 −dF（以及相应力矩）。本处理器把这部分力累加到 TwoPhaseForceBuffers，
 *   由 getTwoPhaseForcesFromLattice() 连同 PSM 动量交换力一起写入 DEM。
 *   原版只传了动量交换力，毛细力对颗粒的作用是缺失的。
 *
 * 【每步调用顺序】（见 oilwaterDrive.cpp 主循环）
 *   setSpheresOnLattice → buffers.prepare() → applyOilWaterShanChen() → collideAndStream
 *   → getTwoPhaseForcesFromLattice() → DEM 子循环
 * ====================================================================================
 */

#ifndef OILWATER_SHAN_CHEN_3D_H
#define OILWATER_SHAN_CHEN_3D_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "plb_ib.h"                 // IBcompositeDynamics / 粒子数据 / 单位换算
#include "liggghtsCouplingWrapper.h"
#include <vector>

namespace plb {

/// 两相力的"中转站"：存颗粒位置（格子单位）和本步累加的毛细反作用力/力矩（格子单位）。
/// 每个 MPI 进程各有一份，只统计本进程格子上算出的贡献；最后与 PSM 力一起归约。
template <typename T>
struct TwoPhaseForceBuffers {
    std::vector<T> capForce;    ///< 3*nPart，毛细反作用力（格子单位）
    std::vector<T> capTorque;   ///< 3*nPart，对应力矩（格子单位）
    typename ParticleData<T>::ParticleDataArrayVector x_lb;   ///< 颗粒中心（格子坐标）

    /// 每个时间步调用一次：更新颗粒位置、把累加数组清零。
    void prepare(LiggghtsCouplingWrapper& wrapper, PhysUnits3D<T> const& units);
};

/// ★★ 2026-10-26（提速）：每步先把两相密度和"固体标记"抽成三个标量场（一次遍历）。
/// 之后 Shan-Chen 处理器读 18 个邻居时直接查这三个场，不再对每个邻居重新累加 19 个分布函数
/// （旧版每个格点要读 18×2×19 = 684 个数，内存访问量约是现在的 20 倍，内核是整步最慢的部分）。
/// 标记：0 = 流体（含固含率 ≤ solidThreshold 的部分覆盖格点），
///       1 = 颗粒（固含率 > solidThreshold，Shan-Chen 把它当固体），2 = 计算域侧壁（BounceBack）。
/// blocks：[0] 水格子 [1] 油格子 [2] ρ_w 场 [3] ρ_o 场 [4] 标记场
template <typename T, template <typename U> class Descriptor>
class ExtractTwoPhaseDensities3D : public BoxProcessingFunctional3D {
public:
    explicit ExtractTwoPhaseDensities3D(T solidThreshold_) : solidThreshold(solidThreshold_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual ExtractTwoPhaseDensities3D<T, Descriptor>* clone() const { return new ExtractTwoPhaseDensities3D<T, Descriptor>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing; modified[1] = modif::nothing;
        modified[2] = modif::staticVariables; modified[3] = modif::staticVariables; modified[4] = modif::staticVariables;
    }
private:
    T solidThreshold;
};

/// 对同一位置的"水格子 + 油格子"施加 Shan-Chen 相互作用、公共速度和润湿性的盒处理器。
/// blocks：[0] 水格子 [1] 油格子 [2] ρ_w 场 [3] ρ_o 场 [4] 标记场（由 ExtractTwoPhaseDensities3D 刚算好）
template <typename T, template <typename U> class Descriptor>
class OilWaterShanChenFunctional3D
    : public BoxProcessingFunctional3D {
public:
    /// G_：跨组分排斥强度；omegaW_/omegaO_：水/油弛豫频率 ω=1/τ；
    /// rhoResidual_：另一相的最低（残余）密度；
    /// wallK_：3 个颗粒类型各自的润湿系数 K = WETTING_GAIN·cotθ（下标 0,1,2 对应 type 1,2,3；θ=90° 时 K=0）；
    /// containerWallK_：计算域四周反弹壁（y/z 侧壁）的 K；
    /// mixTau_/satLow_/satHigh_：局部混合松弛时间开关与插值阈值；uCap_：痕量组分偏离公共速度的上限；
    /// wrapper_：用于查颗粒类型；buffers_：毛细反作用力的累加数组（可为 0 表示不累加）。
    OilWaterShanChenFunctional3D(T G_, T omegaW_, T omegaO_, T rhoResidual_,
                                 Array<T,3> const& wallK_, T containerWallK_,
                                 bool mixTau_, T satLow_, T satHigh_, T uCap_, T tauDiffusion_,
                                 LiggghtsCouplingWrapper* wrapper_,
                                 TwoPhaseForceBuffers<T>* buffers_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);

    virtual OilWaterShanChenFunctional3D<T, Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private:
    /// 二阶平衡态（显式全密度公式，c_s² = 1/3）
    T equilibrium(plint iPop, T rho, Array<T, 3> const& u) const;
    /// 读取一个格点的全密度（Palabos 存的是 f_i − t_i，所以要 +1）和动量
    static void rhoAndJ(Cell<T, Descriptor> const& cell, T& rho, Array<T, 3>& j);
    /// 查某个 LIGGGHTS 颗粒的类型（1/2/3），查不到返回 0
    int particleType(plint partId, plint& liggghtsIndex) const;

    T G, omegaW, omegaO, rhoResidual;
    bool mixTau;          ///< true = 局部混合松弛时间（两相同一 ω），false = 旧做法（各自固定 ω）
    T satLow, satHigh;    ///< 混合粘度插值用的含水比例上下阈值
    T uCap;               ///< 痕量组分"偏离公共速度"的上限（格子单位）
    T omegaD;             ///< 组分间动量交换（扩散）松弛频率 1/τ_D；≤0 表示沿用粘性 ω（旧做法）
    T nuW, nuO;           ///< 水、油格子运动粘度（由 omegaW/omegaO 换算）
    Array<T,3> wallK;     ///< 三类颗粒的润湿系数 K = gain·cotθ
    T containerWallK;     ///< 域侧壁（BounceBack 格点）的 K
    LiggghtsCouplingWrapper* wrapper;
    TwoPhaseForceBuffers<T>* buffers;
};

/// 便捷封装：先抽密度场，再在 domain（一般为流体内部，不含入口面/出口压力层/壁面）上施加两相修正。
/// rhoWField/rhoOField/solidFlagField：由调用者用 MultiScalarField3D<T>(latticeW) 建好、反复使用的工作场。
template <typename T, template <typename U> class Descriptor>
void applyOilWaterShanChen(MultiBlockLattice3D<T, Descriptor>& latticeW,
                           MultiBlockLattice3D<T, Descriptor>& latticeO,
                           Box3D const& domain,
                           T G, T omegaW, T omegaO, T rhoResidual,
                           Array<T,3> const& wallK, T containerWallK,
                           bool mixTau, T satLow, T satHigh, T uCap, T tauDiffusion,
                           LiggghtsCouplingWrapper* wrapper,
                           TwoPhaseForceBuffers<T>* buffers,
                           MultiScalarField3D<T>& rhoWField, MultiScalarField3D<T>& rhoOField,
                           MultiScalarField3D<T>& solidFlagField, T solidThreshold);

}  // namespace plb

#include "oilwaterShanChen3D.hh"

#endif  // OILWATER_SHAN_CHEN_3D_H

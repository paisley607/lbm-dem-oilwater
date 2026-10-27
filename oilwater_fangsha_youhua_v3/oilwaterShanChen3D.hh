/*
 * oilwaterShanChen3D.hh — OilWaterShanChenFunctional3D 与 TwoPhaseForceBuffers 的实现。
 * 物理与公式说明见 oilwaterShanChen3D.h 顶部注释；这里逐步解释代码。
 */

#include "ibDef.h"
#include "utils.h"
#include "lammps.h"
#include "atom.h"
#include <cmath>

namespace plb {

/* ------------------------------------------------------------------------------
 * TwoPhaseForceBuffers：每步开始时刷新颗粒位置并清零累加数组
 * ------------------------------------------------------------------------------ */
template <typename T>
void TwoPhaseForceBuffers<T>::prepare(LiggghtsCouplingWrapper& wrapper,
                                      PhysUnits3D<T> const& units)
{
    // LIGGGHTS 的本地颗粒 + 幽灵颗粒（跨进程边界的副本）都要有位置，
    // 因为本进程格子附近的颗粒可能属于别的进程。
    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    x_lb.resize(nPart);
    for (plint iPart = 0; iPart < nPart; ++iPart) {
        for (int d = 0; d < 3; ++d)
            x_lb[iPart][d] = units.getLbPosition(wrapper.lmp->atom->x[iPart][d]);
    }
    capForce.assign(3*nPart, (T)0);
    capTorque.assign(3*nPart, (T)0);
}

/* ------------------------------------------------------------------------------
 * 构造 / 克隆 / 修改类型声明
 * ------------------------------------------------------------------------------ */
template <typename T, template <typename U> class Descriptor>
OilWaterShanChenFunctional3D<T, Descriptor>::OilWaterShanChenFunctional3D(
        T G_, T omegaW_, T omegaO_, T rhoResidual_,
        Array<T,3> const& wallK_, T containerWallK_,
        bool mixTau_, T satLow_, T satHigh_, T uCap_, T tauDiffusion_,
        LiggghtsCouplingWrapper* wrapper_, TwoPhaseForceBuffers<T>* buffers_)
    : G(G_), omegaW(omegaW_), omegaO(omegaO_), rhoResidual(rhoResidual_),
      mixTau(mixTau_), satLow(satLow_), satHigh(satHigh_), uCap(uCap_),
      omegaD(tauDiffusion_ > (T)0.5 ? (T)1/tauDiffusion_ : (T)0),
      nuW(((T)1/omegaW_ - (T)0.5)/(T)3), nuO(((T)1/omegaO_ - (T)0.5)/(T)3),
      wallK(wallK_), containerWallK(containerWallK_),
      wrapper(wrapper_), buffers(buffers_)
{ }

template <typename T, template <typename U> class Descriptor>
OilWaterShanChenFunctional3D<T, Descriptor>*
OilWaterShanChenFunctional3D<T, Descriptor>::clone() const
{
    return new OilWaterShanChenFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void OilWaterShanChenFunctional3D<T, Descriptor>::getTypeOfModification(
        std::vector<modif::ModifT>& modified) const
{
    // 两套格子的分布函数都被改了，Palabos 需要据此在进程边界同步数据；三个工作场只读。
    modified[0] = modif::staticVariables;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
    modified[4] = modif::nothing;
}

/* ------------------------------------------------------------------------------
 * 抽取两相密度与固体标记（每步一次，供 Shan-Chen 处理器查邻居用）
 * ------------------------------------------------------------------------------ */
template <typename T, template <typename U> class Descriptor>
void ExtractTwoPhaseDensities3D<T, Descriptor>::processGenericBlocks(
        Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    BlockLattice3D<T, Descriptor>& latticeW = *dynamic_cast<BlockLattice3D<T, Descriptor>*>(blocks[0]);
    BlockLattice3D<T, Descriptor>& latticeO = *dynamic_cast<BlockLattice3D<T, Descriptor>*>(blocks[1]);
    ScalarField3D<T>& rhoWF = *dynamic_cast<ScalarField3D<T>*>(blocks[2]);
    ScalarField3D<T>& rhoOF = *dynamic_cast<ScalarField3D<T>*>(blocks[3]);
    ScalarField3D<T>& flagF = *dynamic_cast<ScalarField3D<T>*>(blocks[4]);
    Dot3D const o2 = computeRelativeDisplacement(latticeW, rhoWF);
    Dot3D const o3 = computeRelativeDisplacement(latticeW, rhoOF);
    Dot3D const o4 = computeRelativeDisplacement(latticeW, flagF);
    Dot3D const oO = computeRelativeDisplacement(latticeW, latticeO);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX)
    for (plint iY = domain.y0; iY <= domain.y1; ++iY)
    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
        Cell<T, Descriptor>& cw = latticeW.get(iX, iY, iZ);
        Cell<T, Descriptor> const& co = latticeO.get(iX+oO.x, iY+oO.y, iZ+oO.z);
        T rw = (T)1, ro = (T)1;   // Palabos 存 f_i − t_i，求和后 +1 才是密度
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) { rw += cw[iPop]; ro += co[iPop]; }
        T flag = (T)0;
        if (dynamic_cast<BounceBack<T, Descriptor>*>(&cw.getDynamics())) flag = (T)2;   // 侧壁（BounceBack；用 dynamic_cast 判断，老版 Palabos 也能编译）
        else {
            IBdynamicsParticleData<T, Descriptor>* pd = getParticleDataFromCell<T, Descriptor>(cw);
            if (pd && pd->solidFraction > solidThreshold) flag = (T)1;   // 颗粒（固含率超过阈值，修订七）
        }
        rhoWF.get(iX+o2.x, iY+o2.y, iZ+o2.z) = rw;
        rhoOF.get(iX+o3.x, iY+o3.y, iZ+o3.z) = ro;
        flagF.get(iX+o4.x, iY+o4.y, iZ+o4.z) = flag;
    }
}

/* ------------------------------------------------------------------------------
 * 小工具
 * ------------------------------------------------------------------------------ */
/// ★★ 修订八：痕量组分"密度地板"修复次数（本进程累计；主程序在日志步汇总打印后清零）
inline long& traceDensityRepairCount() { static long n = 0; return n; }
template <typename T, template <typename U> class Descriptor>
T OilWaterShanChenFunctional3D<T, Descriptor>::equilibrium(
        plint iPop, T rho, Array<T, 3> const& u) const
{
    // f_eq_i = t_i ρ [1 + 3 c·u + (9/2)(c·u)² − (3/2) u²]
    T const cDotU = Descriptor<T>::c[iPop][0]*u[0]
                  + Descriptor<T>::c[iPop][1]*u[1]
                  + Descriptor<T>::c[iPop][2]*u[2];
    T const uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    return Descriptor<T>::t[iPop]*rho*((T)1 + (T)3*cDotU + (T)4.5*cDotU*cDotU - (T)1.5*uSqr);
}

template <typename T, template <typename U> class Descriptor>
void OilWaterShanChenFunctional3D<T, Descriptor>::rhoAndJ(
        Cell<T, Descriptor> const& cell, T& rho, Array<T, 3>& j)
{
    // Palabos 为了精度把分布函数存成 f_i − t_i（"去偏移"存储），
    // 直接求和得到的是 ρ−1，所以全密度要 +1；动量不受影响（Σ t_i c_i = 0）。
    T rhoBar = (T)0;
    j[0] = j[1] = j[2] = (T)0;
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rhoBar += cell[iPop];
        for (int d = 0; d < 3; ++d) j[d] += cell[iPop]*Descriptor<T>::c[iPop][d];
    }
    rho = rhoBar + (T)1;
}


template <typename T, template <typename U> class Descriptor>
int OilWaterShanChenFunctional3D<T, Descriptor>::particleType(plint partId, plint& liggghtsIndex) const
{
    liggghtsIndex = -1;
    if (!wrapper || partId < 1) return 0;
    liggghtsIndex = wrapper->lmp->atom->map((int)partId);   // 颗粒编号 → 本进程数组下标
    if (liggghtsIndex < 0) return 0;                          // 本进程没有这个颗粒的副本
    return wrapper->lmp->atom->type[liggghtsIndex];
}

/* ------------------------------------------------------------------------------
 * 核心：逐格点施加两相相互作用
 * ------------------------------------------------------------------------------ */
template <typename T, template <typename U> class Descriptor>
void OilWaterShanChenFunctional3D<T, Descriptor>::processGenericBlocks(
        Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    static const plint q = Descriptor<T>::q;
    BlockLattice3D<T, Descriptor>& latticeW = *dynamic_cast<BlockLattice3D<T, Descriptor>*>(blocks[0]);
    BlockLattice3D<T, Descriptor>& latticeO = *dynamic_cast<BlockLattice3D<T, Descriptor>*>(blocks[1]);
    ScalarField3D<T> const& rhoWF = *dynamic_cast<ScalarField3D<T>*>(blocks[2]);
    ScalarField3D<T> const& rhoOF = *dynamic_cast<ScalarField3D<T>*>(blocks[3]);
    ScalarField3D<T> const& flagF = *dynamic_cast<ScalarField3D<T>*>(blocks[4]);
    Dot3D const oO = computeRelativeDisplacement(latticeW, latticeO);
    Dot3D const o2 = computeRelativeDisplacement(latticeW, rhoWF);
    Dot3D const o3 = computeRelativeDisplacement(latticeW, rhoOF);
    Dot3D const o4 = computeRelativeDisplacement(latticeW, flagF);
    Dot3D const relativePosition = latticeW.getLocation();   // 本块在全局格子里的偏移

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                Cell<T, Descriptor>& cellW = latticeW.get(iX, iY, iZ);
                Cell<T, Descriptor>& cellO = latticeO.get(iX+oO.x, iY+oO.y, iZ+oO.z);

                // ---- 0. 颗粒格点（固含率 > 阈值）和侧壁不施加流体相互作用 ----
                // 它们在 PSM 碰撞里做反弹；它们对周围流体的作用通过"虚拟密度"体现（见下）。
                // ★★ 修订七：颗粒格点仍按局部含水比例设混合粘度（PSM 的 B 加权、部分覆盖格点的流体碰撞都要用 τ）。
                T const flagHere = flagF.get(iX+o4.x, iY+o4.y, iZ+o4.z);
                if (flagHere > (T)1.5) continue;                     // 侧壁
                if (flagHere > (T)0.5) {                             // 颗粒
                    // ★★ 修订八：颗粒覆盖的格点里同样做痕量组分"密度地板"修复（PSM 反弹也可能把痕量组分弄成负值，
                    //   本机测试剖面第 1335 步：骨架颗粒表面固含率 0.82 的格点里水密度 −0.012）
                    {
                        T rwC, roC; Array<T,3> jwC, joC;
                        rhoAndJ(cellW, rwC, jwC); rhoAndJ(cellO, roC, joC);
                        if (rwC < (T)0.5*rhoResidual && roC > rhoResidual) {
                            Array<T,3> uRef; for (int d = 0; d < 3; ++d) uRef[d] = joC[d]/roC;
                            for (plint iPop = 0; iPop < q; ++iPop) cellW[iPop] = equilibrium(iPop, rhoResidual, uRef) - Descriptor<T>::t[iPop];
                            ++traceDensityRepairCount();
                        }
                        if (roC < (T)0.5*rhoResidual && rwC > rhoResidual) {
                            Array<T,3> uRef; for (int d = 0; d < 3; ++d) uRef[d] = jwC[d]/rwC;
                            for (plint iPop = 0; iPop < q; ++iPop) cellO[iPop] = equilibrium(iPop, rhoResidual, uRef) - Descriptor<T>::t[iPop];
                            ++traceDensityRepairCount();
                        }
                    }
                    if (mixTau) {
                        T rwH = rhoWF.get(iX+o2.x, iY+o2.y, iZ+o2.z), roH = rhoOF.get(iX+o3.x, iY+o3.y, iZ+o3.z);
                        if (rwH < rhoResidual) rwH = rhoResidual;
                        if (roH < rhoResidual) roH = rhoResidual;
                        T phiH = (rwH/(rwH + roH) - satLow)/(satHigh - satLow);
                        if (phiH < (T)0) phiH = (T)0;
                        if (phiH > (T)1) phiH = (T)1;
                        T const omH = (T)1/((T)0.5 + (T)3/(phiH/nuW + ((T)1 - phiH)/nuO));
                        T const omTrace = (omegaD > (T)0) ? omegaD : (T)1;          // 痕量组分同上（修订八）
                        T const xWH = rwH/(rwH + roH);
                        T wW = xWH/satLow;            if (wW > (T)1) wW = (T)1;
                        T wO = ((T)1 - xWH)/satLow;   if (wO > (T)1) wO = (T)1;
                        cellW.getDynamics().setOmega(wW*omH + ((T)1 - wW)*omTrace);
                        cellO.getDynamics().setOmega(wO*omH + ((T)1 - wO)*omTrace);
                    }
                    continue;
                }

                // ---- 1. 本格的密度与动量 ----
                T rhoW, rhoO;
                Array<T, 3> jW, jO;
                rhoAndJ(cellW, rhoW, jW);
                rhoAndJ(cellO, rhoO, jO);
                // ★★ 修订八：痕量组分"密度地板"。贴壁、界面附近的痕量组分（油区里的水、水区里的油）偶尔会被
                //   流走到接近 0 甚至负值（本机测试剖面第 1343 步：骨架颗粒 63 表面一排格点的水密度为 −0.001），
                //   之后下面的力计算只能用截断值 ρ_res，给它的动量与真实质量对不上，越来越坏直至发散。
                //   低于 ρ_res/2 时，把这个组分直接重置为"密度 ρ_res、速度与另一组分相同"的平衡态。
                //   每次最多补进 0.005 的痕量质量（约为该格点总质量的 0.5%），对结果的影响可以忽略；
                //   主程序在日志里打印 [REPAIR] 次数，次数持续很多说明参数（G、τ、流速）需要调整。
                if (rhoW < (T)0.5*rhoResidual && rhoO > rhoResidual) {
                    Array<T,3> uRef; for (int d = 0; d < 3; ++d) uRef[d] = jO[d]/rhoO;
                    for (plint iPop = 0; iPop < q; ++iPop) cellW[iPop] = equilibrium(iPop, rhoResidual, uRef) - Descriptor<T>::t[iPop];
                    rhoAndJ(cellW, rhoW, jW);
                    ++traceDensityRepairCount();
                }
                if (rhoO < (T)0.5*rhoResidual && rhoW > rhoResidual) {
                    Array<T,3> uRef; for (int d = 0; d < 3; ++d) uRef[d] = jW[d]/rhoW;
                    for (plint iPop = 0; iPop < q; ++iPop) cellO[iPop] = equilibrium(iPop, rhoResidual, uRef) - Descriptor<T>::t[iPop];
                    rhoAndJ(cellO, rhoO, jO);
                    ++traceDensityRepairCount();
                }
                if (rhoW < rhoResidual) rhoW = rhoResidual;   // 不允许出现 0 密度
                if (rhoO < rhoResidual) rhoO = rhoResidual;

                // ---- 2. 邻居求和 H_σ = Σ_i t_i ρ̃_σ(x+c_i) c_i ----
                // 邻居是流体：用它真实的密度。邻居是固体（颗粒内部格点或计算域侧壁）：先记下方向，
                // 等流体邻居都扫完、算出本格的相场梯度后，再按接触角给它一个"虚拟密度"（见 2b）。
                //
                // ★★ 2026-10-26 修订六：润湿模型改为"几何接触角"（Ding & Spelt 2007, PRE 75:046708；
                //   用于伪势/Shan-Chen 模型见 Li, Luo, Kang, Chen 2014, PRE 90:053301）。
                //   旧的"虚拟密度"做法（修订三：固体对非润湿相施加 (2s−1)ρ_wall 的排斥）在整个固体表面、
                //   包括远离油水界面的本体流体上都施力：非润湿相被从壁面推开、在窄喉道中间被挤压，
                //   numpy 二维多孔介质算例里 s=0.3（油湿）和 s≥0.6（水湿）都在 500 步内发散。
                //   几何接触角只在"三相接触线"附近起作用：
                //     φ = (ρ_w−ρ_o)/(ρ_w+ρ_o) 为相场，n 为由固体指向流体的单位法向，
                //     壁面处要求 n·∇φ = −cotθ · |∇φ − (n·∇φ)n|（θ 为从水相一侧量的接触角）。
                //   据此把固体邻格（位于 x+c_i）的 φ 外推为
                //     φ_s = φ(x) + c_i·∇_tφ − K·|∇_tφ|·(c_i·n)，K = WETTING_GAIN·cotθ，
                //   总密度取本格的 ρ_w+ρ_o（"镜像"），再按 φ_s 分给两相。本体流体里 ∇φ≈0，φ_s=φ(x)，
                //   壁面对流体没有任何附加力——不再有贴壁的密度亏损层，窄喉道不再被挤爆；
                //   只有接触线附近按 θ 把界面"掰"成需要的角度。
                Array<T, 3> Hw((T)0, (T)0, (T)0), Ho((T)0, (T)0, (T)0);
                Array<T, 3> gradAcc((T)0, (T)0, (T)0);   // Σ t_i φ(x+c_i) c_i（只含流体邻格）
                Array<T, 3> S((T)0, (T)0, (T)0);         // Σ t_i c_i（只含固体邻格），−S 即壁面法向
                T const rhoT = rhoW + rhoO;
                T const phiHere = (rhoW - rhoO)/rhoT;
                int nSolid = 0;
                plint solidPop[18];      // 固体邻格的方向
                T solidK[18];            // 该方向固体的 K = gain·cotθ
                plint solidInd[18];      // 颗粒在 LIGGGHTS 本地数组的下标（侧壁为 −1）

                for (plint iPop = 1; iPop < q; ++iPop) {
                    plint const nX = iX + Descriptor<T>::c[iPop][0];
                    plint const nY = iY + Descriptor<T>::c[iPop][1];
                    plint const nZ = iZ + Descriptor<T>::c[iPop][2];
                    T const tw = Descriptor<T>::t[iPop];

                    // ★★ 修订四：计算域侧壁（y/z 四个面，BounceBack 动力学）是固体，不能读它存的分布函数
                    //   （BounceBack 格点里的"密度"没有物理意义，旧版据此把侧壁当成油、水密度都为 1 的流体，
                    //   第 7 步在侧壁旁发散）。
                    T const flagN = flagF.get(nX+o4.x, nY+o4.y, nZ+o4.z);
                    bool const isWall = flagN > (T)1.5;
                    bool const isParticle = !isWall && flagN > (T)0.5;

                    if (isWall || isParticle) {
                        plint ind = -1;
                        T K = containerWallK;
                        if (isParticle) {
                            IBdynamicsParticleData<T, Descriptor>* pdN =
                                getParticleDataFromCell<T, Descriptor>(latticeW.get(nX, nY, nZ));
                            int const type = pdN ? particleType(pdN->partId, ind) : 0;
                            K = (type >= 1 && type <= 3) ? wallK[type-1] : wallK[0];   // 查不到类型按 type 1
                        }
                        solidPop[nSolid] = iPop; solidK[nSolid] = K; solidInd[nSolid] = ind; ++nSolid;
                        for (int d = 0; d < 3; ++d) S[d] += tw*(T)Descriptor<T>::c[iPop][d];
                    } else {
                        T rwN = rhoWF.get(nX+o2.x, nY+o2.y, nZ+o2.z);
                        T roN = rhoOF.get(nX+o3.x, nY+o3.y, nZ+o3.z);
                        if (rwN < (T)0) rwN = (T)0;
                        if (roN < (T)0) roN = (T)0;
                        T const rtN = rwN + roN;
                        T const phiN = (rtN > (T)1e-12) ? (rwN - roN)/rtN : phiHere;
                        for (int d = 0; d < 3; ++d) {
                            T const cc = (T)Descriptor<T>::c[iPop][d];
                            Hw[d] += tw*rwN*cc;
                            Ho[d] += tw*roN*cc;
                            gradAcc[d] += tw*phiN*cc;
                        }
                    }
                }

                // ---- 2b. 固体邻格的虚拟密度（几何接触角）与毛细反作用力 ----
                if (nSolid > 0) {
                    // 相场梯度 ∇φ ≈ 3 Σ t_i φ(x+c_i) c_i，固体方向暂按 φ(x)（镜像）计
                    Array<T,3> grad, nrm;
                    for (int d = 0; d < 3; ++d) grad[d] = (T)3*(gradAcc[d] + phiHere*S[d]);
                    T const Sn = std::sqrt(S[0]*S[0] + S[1]*S[1] + S[2]*S[2]);
                    // 法向：由固体指向流体。两侧都是固体（1 格宽缝隙）时 S≈0，法向无定义 → 只做镜像
                    for (int d = 0; d < 3; ++d) nrm[d] = (Sn > (T)1e-12) ? -S[d]/Sn : (T)0;
                    T const gn = grad[0]*nrm[0] + grad[1]*nrm[1] + grad[2]*nrm[2];
                    Array<T,3> gt;
                    for (int d = 0; d < 3; ++d) gt[d] = grad[d] - gn*nrm[d];
                    T const gtm = std::sqrt(gt[0]*gt[0] + gt[1]*gt[1] + gt[2]*gt[2]);
                    T const phiMax = (T)1 - (T)2*rhoResidual;
                    T const xG = (T)(relativePosition.x + iX);
                    T const yG = (T)(relativePosition.y + iY);
                    T const zG = (T)(relativePosition.z + iZ);
                    for (int k = 0; k < nSolid; ++k) {
                        plint const iPop = solidPop[k];
                        T const tw = Descriptor<T>::t[iPop];
                        T cgt = 0, cn = 0;
                        for (int d = 0; d < 3; ++d) {
                            cgt += (T)Descriptor<T>::c[iPop][d]*gt[d];
                            cn  += (T)Descriptor<T>::c[iPop][d]*nrm[d];
                        }
                        T phiS = phiHere + cgt - solidK[k]*gtm*cn;    // c_i 指向固体，cn<0：水湿(K>0) 时 φ_s 偏向水
                        if (phiS >  phiMax) phiS =  phiMax;
                        if (phiS < -phiMax) phiS = -phiMax;
                        T const rwS = (T)0.5*rhoT*((T)1 + phiS);      // 水"在壁面上"的虚拟密度
                        T const roS = (T)0.5*rhoT*((T)1 - phiS);      // 油"在壁面上"的虚拟密度
                        for (int d = 0; d < 3; ++d) {
                            T const cc = (T)Descriptor<T>::c[iPop][d];
                            Hw[d] += tw*rwS*cc;
                            Ho[d] += tw*roS*cc;
                        }
                        // 毛细反作用力：流体在本格受固体方向的力 dF = −G(ρ_w ρ̃_o + ρ_o ρ̃_w) t_i c_i，
                        // 颗粒受 −dF；力矩以"本格 − 颗粒中心"为力臂（格子单位）。侧壁（ind<0）不记。
                        plint const ind = solidInd[k];
                        if (buffers && ind >= 0 && (plint)buffers->x_lb.size() > ind) {
                            T const mag = G*(rhoW*roS + rhoO*rwS)*tw;
                            Array<T,3> dF;
                            for (int d = 0; d < 3; ++d) dF[d] = mag*(T)Descriptor<T>::c[iPop][d];
                            T const lx = xG - buffers->x_lb[ind][0];
                            T const ly = yG - buffers->x_lb[ind][1];
                            T const lz = zG - buffers->x_lb[ind][2];
                            buffers->capForce[3*ind+0] += dF[0];
                            buffers->capForce[3*ind+1] += dF[1];
                            buffers->capForce[3*ind+2] += dF[2];
                            buffers->capTorque[3*ind+0] += ly*dF[2] - lz*dF[1];
                            buffers->capTorque[3*ind+1] += lz*dF[0] - lx*dF[2];
                            buffers->capTorque[3*ind+2] += lx*dF[1] - ly*dF[0];
                        }
                    }
                }

                // ---- 3. 跨组分力 ----
                Array<T, 3> Fw, Fo;
                for (int d = 0; d < 3; ++d) {
                    Fw[d] = -G*rhoW*Ho[d];
                    Fo[d] = -G*rhoO*Hw[d];
                }

                // ---- 4. 公共速度 u'（两相粘性耦合的来源）----
                // ---- 4a. 本格两相用的松弛频率 ----
                // mixTau：按局部含水比例插值混合粘度，两套格子用同一个 ω（并写回格点动力学，
                //         这样随后 collideAndStream 的碰撞和 PSM 的 B 加权都用它）。
                // 否则：沿用水/油各自固定的 ω（旧做法，只在 M≲9 时稳定）。
                T omW = omegaW, omO = omegaO;
                if (mixTau) {
                    T const sat = rhoW/(rhoW + rhoO);
                    T phi = (sat - satLow)/(satHigh - satLow);
                    if (phi < (T)0) phi = (T)0;
                    if (phi > (T)1) phi = (T)1;
                    T const nuMix = (T)1/(phi/nuW + ((T)1 - phi)/nuO);
                    T const omMix = (T)1/((T)0.5 + (T)3*nuMix);
                    omW = omO = omMix;
                    // ★★ 修订八：痕量组分的应力用 τ=τ_D（默认 1）松弛，主组分仍用混合 τ。
                    //   油区 τ≈5.45 时，痕量水的非平衡应力衰减很慢，而它的速度每步都被动量交换"拉"向公共速度，
                    //   非平衡部分相对它那一点点密度会变得很大，分布函数变负、密度跟着变负（第 1343 步发散的起因）。
                    //   痕量组分只占该格点质量的不到 1%，它的粘度对混合物粘度没有可见影响。
                    //   做法：含量低于 MIX_SAT_LOW（5%）时按含量线性过渡到 τ_D，界面区（两相都 >5%）不受影响。
                    T const omTrace = (omegaD > (T)0) ? omegaD : (T)1;
                    T const xW = rhoW/(rhoW + rhoO);
                    T wW = xW/satLow;            if (wW > (T)1) wW = (T)1;
                    T wO = ((T)1 - xW)/satLow;   if (wO > (T)1) wO = (T)1;
                    cellW.getDynamics().setOmega(wW*omMix + ((T)1 - wW)*omTrace);
                    cellO.getDynamics().setOmega(wO*omMix + ((T)1 - wO)*omTrace);
                }
                // ---- 4b. 组分间动量交换用的松弛频率（修订五）----
                // omegaD>0：两相都用 ω_D=1/τ_D（默认 τ_D=1），u' 就是质量加权的公共速度，两相动量严格守恒；
                // 粘性仍由上面的 ω_mix 在随后的正则化碰撞里体现（setOmega 已写回）。
                // omegaD=0（DIFFUSION_TAU=0）：沿用粘性 ω（旧做法）。
                if (omegaD > (T)0) { omW = omegaD; omO = omegaD; }
                T const denom = omW*rhoW + omO*rhoO;
                Array<T, 3> uPrime;
                for (int d = 0; d < 3; ++d) uPrime[d] = (omW*jW[d] + omO*jO[d])/denom;

                // ---- 5. 目标动量 j* = (1−ω) j + ω ρ u' + F ----
                // 正则化碰撞会保持动量，所以碰撞后的动量就是 j*，每步净注入恰好为 F。
                Array<T, 3> uOldW, uNewW, uOldO, uNewO;
                for (int d = 0; d < 3; ++d) {
                    T const jWstar = ((T)1 - omW)*jW[d] + omW*rhoW*uPrime[d] + Fw[d];
                    T const jOstar = ((T)1 - omO)*jO[d] + omO*rhoO*uPrime[d] + Fo[d];
                    uOldW[d] = jW[d]/rhoW;   uNewW[d] = jWstar/rhoW;
                    uOldO[d] = jO[d]/rhoO;   uNewO[d] = jOstar/rhoO;
                }
                // ---- 5'. 痕量组分"偏离公共速度"的上限（2026-10-25 修订二）----
                // 界面两侧、壁面旁只有残余密度的那个组分受到单侧排斥，F/ρ 与它自己的密度无关，
                // 目标速度可能超过 1 而毫无物理意义，却会让分布函数变负、计算发散，所以要限幅。
                // 旧版把整个目标速度缩放 → 连切向（沿界面流动方向）的速度也被缩小，
                // 界面每步损失动量、变成"无滑移墙"（M=1 分层流误差 96%）。
                // 现在（mixTau=true 时）只缩放 v = u_new − u' 这一"偏离公共速度"的部分，u' 本身原样保留。
                // 主相的偏离量 ~1e-4，完全不受影响。
                // mixTau=false（两套格子各用固定 τ）时保留旧的整体缩放：那种做法本身在壁面附近
                // 不够稳定，旧的整体缩放带来的额外阻尼恰好是它在 τ_o≲3.5 时能算下去的原因（numpy 验证：
                // 固定 τ + 只限偏离量，在 s_wall≥0.5 的接触角算例里发散）。
                if (mixTau) {
                    Array<T,3> vW, vO;
                    for (int d = 0; d < 3; ++d) { vW[d] = uNewW[d] - uPrime[d]; vO[d] = uNewO[d] - uPrime[d]; }
                    T const sW = std::sqrt(vW[0]*vW[0] + vW[1]*vW[1] + vW[2]*vW[2]);
                    T const sO = std::sqrt(vO[0]*vO[0] + vO[1]*vO[1] + vO[2]*vO[2]);
                    if (sW > uCap) for (int d = 0; d < 3; ++d) uNewW[d] = uPrime[d] + vW[d]*uCap/sW;
                    if (sO > uCap) for (int d = 0; d < 3; ++d) uNewO[d] = uPrime[d] + vO[d]*uCap/sO;
                } else {
                    T const sW = std::sqrt(uNewW[0]*uNewW[0] + uNewW[1]*uNewW[1] + uNewW[2]*uNewW[2]);
                    T const sO = std::sqrt(uNewO[0]*uNewO[0] + uNewO[1]*uNewO[1] + uNewO[2]*uNewO[2]);
                    if (sW > uCap) for (int d = 0; d < 3; ++d) uNewW[d] *= uCap/sW;
                    if (sO > uCap) for (int d = 0; d < 3; ++d) uNewO[d] *= uCap/sO;
                }

                // ---- 6. 只改平衡部分、保留非平衡部分：f ← f + feq(ρ,u_new) − feq(ρ,u_old) ----
                for (plint iPop = 0; iPop < q; ++iPop) {
                    cellW[iPop] += equilibrium(iPop, rhoW, uNewW) - equilibrium(iPop, rhoW, uOldW);
                    cellO[iPop] += equilibrium(iPop, rhoO, uNewO) - equilibrium(iPop, rhoO, uOldO);
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------
 * 便捷封装
 * ------------------------------------------------------------------------------ */
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
                           MultiScalarField3D<T>& solidFlagField, T solidThreshold)
{
    std::vector<MultiBlock3D*> blocks;
    blocks.push_back(&latticeW);
    blocks.push_back(&latticeO);
    blocks.push_back(&rhoWField);
    blocks.push_back(&rhoOField);
    blocks.push_back(&solidFlagField);
    // (a) 全场抽两相密度和固体标记（Palabos 会自动同步各进程交界处的 envelope）
    applyProcessingFunctional(new ExtractTwoPhaseDensities3D<T, Descriptor>(solidThreshold), latticeW.getBoundingBox(), blocks);
    // (b) 在流体内部施加两相相互作用
    applyProcessingFunctional(
        new OilWaterShanChenFunctional3D<T, Descriptor>(G, omegaW, omegaO, rhoResidual,
                                                        wallK, containerWallK,
                                                        mixTau, satLow, satHigh, uCap, tauDiffusion,
                                                        wrapper, buffers),
        domain, blocks);
}

}  // namespace plb

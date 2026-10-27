/*
 * diagScan3D.h —— 排错用的"异常格点扫描"（2026-10-26 新增）
 * ====================================================================================
 * 【用途】计算一开始就发散时，想知道"最先出问题的是哪里"：侧壁旁、入口、出口、颗粒表面、
 *   颗粒内部还是孔隙中间。参数 DIAGNOSTIC_SCAN_STEPS=N（默认 0 = 关闭）时，前 N 步每步碰撞+迁移后
 *   扫描全场，把可疑格点（总密度偏离 1 超过 12%、出现负密度、混合物速度 > 0.02 格子单位、NaN）
 *   按位置分类统计个数，并打印每类最严重的 4 个格点（每进程每步），形如：
 *     [DIAG] step=0 cls=wallAdj x=4 y=1 z=3 sf=0.000 pid=0 rw=0.8892 ro=0.2109 rt=1.1001 |u|=0.15 ...
 *   cls 分类：bulk=孔隙中间  wallAdj=离侧壁 ≤2 格  inlet=x≤3  outlet=x≥nx−4
 *             psmSurf=颗粒表面（部分覆盖）  psmInner=颗粒内部  bbWall=侧壁格点本身
 *   sf=固含率，pid=颗粒编号，rw/ro/rt=水/油/总密度，|u|=混合物速度（格子单位）。
 * 【注意】只读不写，不影响计算结果；但扫描全场要花时间，正常计算请保持 0。
 * ====================================================================================
 */
#ifndef OILWATER_DIAG_SCAN_3D_H
#define OILWATER_DIAG_SCAN_3D_H

#include <cstdio>
#include <cmath>

namespace plb {

/// 每步、每类"最可疑"的前 4 个格点（每个 MPI 进程一份），扫描完再统一打印
struct DiagWorst {
    int count[8]; double score[8][4]; char line[8][4][320];
    DiagWorst() { for (int i=0;i<8;++i) { count[i]=0; for (int k=0;k<4;++k) { score[i][k]=-1.; line[i][k][0]='\0'; } } }
    void offer(int cls, double sc, const char* txt) {
        ++count[cls];
        int k=3; if (sc<=score[cls][k]) return;
        while (k>0 && sc>score[cls][k-1]) { score[cls][k]=score[cls][k-1]; std::snprintf(line[cls][k],320,"%s",line[cls][k-1]); --k; }
        score[cls][k]=sc; std::snprintf(line[cls][k],320,"%s",txt);
    }
};
inline DiagWorst& diagWorst() { static DiagWorst w; return w; }
static const char* const diagClassName[8] = {"bulk","wallAdj","inlet","outlet","psmSurf","psmInner","bbWall","other"};

template<typename T, template<typename U> class Descriptor>
class DiagScan3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor> {
public:
    DiagScan3D(plint step_, plint nx_, plint ny_, plint nz_, int maxPer_)
        : step(step_), nx(nx_), ny(ny_), nz(nz_), maxPer(maxPer_) { }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lw, BlockLattice3D<T,Descriptor>& lo) {
        Dot3D const off = lw.getLocation();
        for (plint x=domain.x0; x<=domain.x1; ++x)
        for (plint y=domain.y0; y<=domain.y1; ++y)
        for (plint z=domain.z0; z<=domain.z1; ++z) {
            Cell<T,Descriptor>& cw = lw.get(x,y,z);
            Cell<T,Descriptor>& co = lo.get(x,y,z);
            // 全密度与动量（Palabos 存 f−t，求和后 +1）
            T rw=1, ro=1; Array<T,3> jw, jo; jw.resetToZero(); jo.resetToZero();
            for (plint i=0;i<Descriptor<T>::q;++i) {
                rw += cw[i]; ro += co[i];
                for (int d=0;d<3;++d) { jw[d]+=cw[i]*Descriptor<T>::c[i][d]; jo[d]+=co[i]*Descriptor<T>::c[i][d]; }
            }
            plint const gx=x+off.x, gy=y+off.y, gz=z+off.z;
            T sf=0; plint pid=0;
            IBdynamicsParticleData<T,Descriptor>* pd = getParticleDataFromCell<T,Descriptor>(cw);
            if (pd) { sf=pd->solidFraction; pid=pd->partId; }
            int cls=0;
            if (dynamic_cast<BounceBack<T,Descriptor>*>(&cw.getDynamics())) cls=6;
            else if (sf>(T)0.999) cls=5;
            else if (sf>(T)0) cls=4;
            else if (gx<=3) cls=2;
            else if (gx>=nx-4) cls=3;
            else if (gy<=2||gy>=ny-3||gz<=2||gz>=nz-3) cls=1;
            T const rt=rw+ro;
            T um[3]; for (int d=0;d<3;++d) um[d]=(jw[d]+jo[d])/rt;
            T const us=std::sqrt(um[0]*um[0]+um[1]*um[1]+um[2]*um[2]);
            bool const bad = !(rt==rt) || !(us==us) || rt>(T)1.12 || rt<(T)0.88 || rw<(T)0 || ro<(T)0 || us>(T)0.02;
            if (bad) {
                // 可疑程度：NaN 最高；其次负密度；再按 |ρ_总−1|/0.12 与 |u|/0.02 的较大者
                double sc = (!(rt==rt) || !(us==us)) ? 1e30 : std::max(std::fabs((double)rt-1.)/0.12, (double)us/0.02);
                if (rw<(T)0 || ro<(T)0) sc += 1e6;
                char txt[320];
                std::snprintf(txt,320,"[DIAG] step=%ld cls=%s x=%ld y=%ld z=%ld sf=%.3f pid=%ld rw=%.4f ro=%.4f rt=%.4f |u|=%.4g u=(%.3g,%.3g,%.3g)",
                    (long)step, diagClassName[cls], (long)gx,(long)gy,(long)gz,(double)sf,(long)pid,
                    (double)rw,(double)ro,(double)rt,(double)us,(double)um[0],(double)um[1],(double)um[2]);
                diagWorst().offer(cls, sc, txt);
            }
        }
    }
    virtual DiagScan3D* clone() const { return new DiagScan3D(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& m) const { m[0]=modif::nothing; m[1]=modif::nothing; }
private:
    plint step, nx, ny, nz; int maxPer;
};

/// 主循环里调用：清零计数后扫描全场
template<typename T, template<typename U> class Descriptor>
void diagnosticScan(MultiBlockLattice3D<T,Descriptor>& latticeW, MultiBlockLattice3D<T,Descriptor>& latticeO, plint step)
{
    diagWorst() = DiagWorst();
    std::vector<MultiBlock3D*> a; a.push_back(&latticeW); a.push_back(&latticeO);
    Box3D const bb = latticeW.getBoundingBox();
    applyProcessingFunctional(new DiagScan3D<T,Descriptor>(step, bb.getNx(), bb.getNy(), bb.getNz(), 4), bb, a);
    DiagWorst const& w = diagWorst();
    for (int c=0;c<8;++c) {
        if (w.count[c]==0) continue;
        std::printf("[DIAG] step=%ld cls=%s 可疑格点数(本进程)=%d，最严重的几个：\n", (long)step, diagClassName[c], w.count[c]);
        for (int k=0;k<4;++k) if (w.score[c][k]>=0.) std::printf("%s\n", w.line[c][k]);
    }
    std::fflush(stdout);
}

}  // namespace plb

#endif  // OILWATER_DIAG_SCAN_3D_H

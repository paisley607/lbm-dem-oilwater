/*
 * oilwaterCoupling3D.hh — getTwoPhaseForcesFromLattice 的实现。
 */

#include "lammps.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

#include <vector>
#include <algorithm>
#include <cmath>

namespace plb {

template <typename T, template <typename U> class Descriptor>
plint getTwoPhaseForcesFromLattice(MultiBlockLattice3D<T, Descriptor>& latticeW,
                                   MultiBlockLattice3D<T, Descriptor>& latticeO,
                                   LiggghtsCouplingWrapper& wrapper,
                                   PhysUnits3D<T> const& units,
                                   TwoPhaseForceBuffers<T> const* capillary,
                                   T& maxCapillaryForcePhys)
{
    // static：数组在整个运行期间复用，避免每步重新分配内存
    static std::vector<T> force, torque;
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;

    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    plint const n_force = nPart*3;
    maxCapillaryForcePhys = (T)0;
    if (nPart == 0) return 0;

    // ---- 1. 颗粒位置（格子单位）：SumForceTorque3D 用它算力矩的力臂 ----
    x_lb.resize(nPart);
    for (plint iPart = 0; iPart < nPart; ++iPart)
        for (int d = 0; d < 3; ++d)
            x_lb[iPart][d] = units.getLbPosition(wrapper.lmp->atom->x[iPart][d]);

    // ---- 2. 力/力矩累加数组清零（只清一次）----
    force.assign(n_force, (T)0);
    torque.assign(n_force, (T)0);

    // 本进程恰好一个格子块（LatticeDecomposition 保证）
    plint const iBlock = latticeW.getLocalInfo().getBlocks()[0];
    std::map<plint, Box3D> blockmap = latticeW.getSparseBlockStructure().getBulks();
    Box3D const localBB = blockmap[iBlock];

    // ---- 3. 水格子 PSM 力：逐格把 particleData.hydrodynamicForce 累加到所属颗粒 ----
    applyProcessingFunctional(
        new SumForceTorque3D<T, Descriptor>(x_lb, &force.front(), &torque.front(), wrapper),
        localBB, latticeW);
    // ---- 4. 油格子 PSM 力：继续累加到同一数组 ----
    applyProcessingFunctional(
        new SumForceTorque3D<T, Descriptor>(x_lb, &force.front(), &torque.front(), wrapper),
        localBB, latticeO);

    // ---- 5. 毛细反作用力（格子单位），由 Shan-Chen 处理器累加好 ----
    if (capillary && (plint)capillary->capForce.size() == n_force) {
        for (plint i = 0; i < n_force; ++i) {
            force[i]  += capillary->capForce[i];
            torque[i] += capillary->capTorque[i];
        }
        for (plint iPart = 0; iPart < nPart; ++iPart) {
            T const fx = capillary->capForce[3*iPart+0], fy = capillary->capForce[3*iPart+1],
                    fz = capillary->capForce[3*iPart+2];
            maxCapillaryForcePhys = std::max(maxCapillaryForcePhys,
                                             units.getPhysForce(std::sqrt(fx*fx + fy*fy + fz*fz)));
        }
    }

    // ---- 6. 取耦合 fix，清零一次 ----
    LAMMPS_NS::FixLbCouplingOnetoone* couplingFix =
        dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>(
            wrapper.lmp->modify->find_fix_style("couple/lb/onetoone", 0));
    double** f_liggghts = couplingFix->get_force_ptr();
    double** t_liggghts = couplingFix->get_torque_ptr();
    for (plint iPart = 0; iPart < nPart; ++iPart)
        for (int j = 0; j < 3; ++j) { f_liggghts[iPart][j] = 0; t_liggghts[iPart][j] = 0; }

    // ---- 7. 转物理单位，按颗粒编号写入 ----
    for (plint iPart = 0; iPart < nPart; ++iPart) {
        int const tag = wrapper.lmp->atom->tag[iPart];
        int const liggghts_ind = wrapper.lmp->atom->map(tag);
        for (int j = 0; j < 3; ++j) {
            f_liggghts[liggghts_ind][j] += units.getPhysForce(force[3*iPart+j]);
            t_liggghts[liggghts_ind][j] += units.getPhysTorque(torque[3*iPart+j]);
        }
    }
    couplingFix->comm_force_torque();   // 跨进程把同一颗粒的力加起来（一次 MPI 归约）
    return nPart;
}

}  // namespace plb

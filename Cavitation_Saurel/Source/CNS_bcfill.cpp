
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

// HHN21: vibrating bc
#include <math.h>

using namespace amrex;

struct CnsFillExtDir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
            // HHN21: vibrating bc 
            // ONLY need to implement bc on the high side of d=0 / x-direction
            // bc should be reflective (inversed sign) and vibrating
            for(unsigned int d = 0; d < amrex::SpaceDim; ++d) // main for loop
            {
            if(bcr->hi(d) == amrex::BCType::ext_dir && iv[d] > geom.Domain().bigEnd(d))
            {
            // Compute a dimensionally independent location vector
            amrex::GpuArray<int, 3> loc = {iv[0],
            (amrex::SpaceDim < 2 ? 0 : iv[amrex::min(1, amrex::SpaceDim-1)]),
            (amrex::SpaceDim < 3 ? 0 : iv[amrex::min(2, amrex::SpaceDim-1)])}; // IF spacedim is 1, j = 0 & k = 0

            // In 1D, the boundary must be grounded everywhere, if anywhere
            dest(loc[0], loc[1], loc[2], dcomp) = - 2 * (dest(loc[0], loc[1], loc[2], 1) + dest(loc[0], loc[1], loc[2], 2)) * ( M_PI * 2e6 * 5e-6 * cos(2 * M_PI * 2e6 * time) * (0.5 + 0.5 * cos(M_PI * 2e-3))) - (dest(loc[0], loc[1], loc[2], 1) + dest(loc[0], loc[1], loc[2], 2)) * dest(loc[0], loc[1], loc[2], 11);
            }
            }
        }
};

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the descriptor set up in CNS::variableSetUp.

void cns_bcfill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
    // std::cout<<time<<std::endl;
    GpuBndryFuncFab<CnsFillExtDir> gpu_bndry_func(CnsFillExtDir{});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}

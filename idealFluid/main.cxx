#include "idealFluid.hxx"
#define CCTK_ARGUMENTS
#define CCTK_PASS_CTOC

using namespace AsterX;

extern template void Con2Prim_2DNRNoble<idealFluid>(
    CCTK_INT max_iter, CCTK_REAL tolf,
    idealFluid &plasma);

/***************************************************************************
 * AsterX_Con2Prim
 * ---
 *  Routines implemented:
 *   1) 2DNRNoble
 *
 *   Based on con2primFactory (https://github.com/fedelopezar/con2primFactory)
 *   ****************************************************************************/

template <typename typeEoS>
void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS)
{
    CCTK_INT max_iter = 100;
    CCTK_REAL c2p_tol = 1E-10;
    CCTK_REAL gamma = 4. / 3.;

    CCTK_REAL g_up[4][4];
    CCTK_REAL g_lo[4][4];
    CCTK_REAL cons[NCONS];

    cons[D] = 0.0;
    cons[S1_COV] = 0.0;
    cons[S2_COV] = 0.0;
    cons[S3_COV] = 0.0;
    cons[TAU] = 0.0;
    cons[B1] = 0.0;
    cons[B2] = 0.0;
    cons[B3] = 0.0;

    CCTK_REAL prims[NPRIMS];
    prims[RHO] = 0.0;
    prims[V1_CON] = 0.0;
    prims[V2_CON] = 0.0;
    prims[V3_CON] = 0.0;
    prims[EPS] = 0.0;

    // Construct con2primFactory object:
    typeEoS plasma_0(gamma, cons, prims, g_lo, g_up);
    // 1) Try 2DNRNoble
    Con2Prim_2DNRNoble(max_iter, c2p_tol, plasma_0);
}

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS)
{
    if (1)
    { // Use this if for idealFluid/tabeos
        AsterX_Con2Prim_typeEoS<idealFluid>(CCTK_PASS_CTOC);
        // CCTK_PASS_CTOC == cctkGH, and more. Preferred over just cctkGH.
    }
}

int main(int argn, char **argv)
{

    return 0;
}

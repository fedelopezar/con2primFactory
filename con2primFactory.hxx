#ifndef CON2PRIM_HXX
#define CON2PRIM_HXX
//#include "AMReX_GpuQualifiers.H"
#define CCTK_DEVICE
//#define CCTK_DEVICE AMREX_GPU_DEVICE
#define CCTK_HOST
//#define CCTK_HOST AMREX_GPU_HOST
//#include "utils.hxx"
#define CCTK_REAL double
#define CCTK_INT int

#include <math.h>

namespace AsterX
{

    constexpr int NCONS = 8;
    constexpr int NPRIMS = 5;

    /* Abstract class con2primFactory */
    class con2primFactory
    {
    public:
        /* The constructor must initialize the following vectors */
        // std::vector<PrimitiveLabel> PrimitiveLabels;
        CCTK_REAL ConservedVars[NCONS];      // Conserved to solve
        CCTK_REAL PrimitiveVarsSeed[NPRIMS]; // Primitive seeds
        CCTK_REAL PrimitiveVars[NPRIMS];     // Primitive solution

        /* These must be set for 2DNRNoble scheme */
        CCTK_INT Failed_2DNRNoble;
        CCTK_INT Nit_2DNRNoble;
        CCTK_REAL W_Seed, vsq_Sol, Ssq, Press_Seed, Z_Seed, Z_Sol, vsq_Seed, bsq;
        CCTK_REAL gcov[10], gcon[10];
        CCTK_REAL Bsq;
        CCTK_REAL BiSi;
        CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();  // From cons (exact)
        CCTK_HOST CCTK_DEVICE void get_Press_Seed(); // From seed prims and cons
        CCTK_HOST CCTK_DEVICE void get_Z_Seed();     // From seed prims and cons
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE void WZ2Prim();
    };

    class idealFluid : public con2primFactory
    {
    public:
        /* Some attributes */
        CCTK_REAL GammaIdealFluid;
        /* Constructor */
        CCTK_HOST CCTK_DEVICE idealFluid(CCTK_REAL gamma, CCTK_REAL (&cons)[NCONS], CCTK_REAL (&prim)[NPRIMS], CCTK_REAL (&gcov)[4][4], CCTK_REAL (&gcon)[4][4], CCTK_INT forTesting);
        CCTK_HOST CCTK_DEVICE void get_ConservedVarsFromPrimVector(CCTK_REAL (&prim)[NPRIMS]);

        /* Called by 2DNRNoble */
        CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();
        CCTK_HOST CCTK_DEVICE void get_Press_Seed();
        CCTK_HOST CCTK_DEVICE void get_Z_Seed();
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
        CCTK_HOST CCTK_DEVICE void WZ2Prim();

        /* Destructor */
        CCTK_HOST CCTK_DEVICE ~idealFluid();
    };
}

#endif // #ifndef CON2PRIM_HXX

template <typename typeEoS>
CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble(
    CCTK_INT max_iter, CCTK_REAL tolf,
    typeEoS &plasma);
#include "con2primFactory.hxx"

namespace AsterX
{

    constexpr int D = 0;
    constexpr int S1_COV = 1;
    constexpr int S2_COV = 2;
    constexpr int S3_COV = 3;
    constexpr int TAU = 4;
    constexpr int B1 = 5;
    constexpr int B2 = 6;
    constexpr int B3 = 7;

    constexpr int RHO = 0;
    constexpr int V1_CON = 1;
    constexpr int V2_CON = 2;
    constexpr int V3_CON = 3;
    constexpr int V1_COV = 1;
    constexpr int V2_COV = 2;
    constexpr int V3_COV = 3;
    constexpr int EPS = 4;

    constexpr int TT = 0;
    constexpr int TX = 1;
    constexpr int TY = 2;
    constexpr int TZ = 3;
    constexpr int XX = 4;
    constexpr int XY = 5;
    constexpr int XZ = 6;
    constexpr int YY = 7;
    constexpr int YZ = 8;
    constexpr int ZZ = 9;

    // class idealFluid : public con2primFactory
    // {
    // public:
    //     /* Some attributes */
    //     CCTK_REAL GammaIdealFluid;
    //     /* Constructor */
    //     CCTK_HOST CCTK_DEVICE idealFluid(CCTK_REAL gamma, CCTK_REAL (&cons)[NCONS], CCTK_REAL (&prim)[NPRIMS], CCTK_REAL (&gcov)[4][4], CCTK_REAL (&gcon)[4][4]);

    //     /* Called by 2DNRNoble */
    //     CCTK_HOST CCTK_DEVICE void get_Ssq_Exact();
    //     CCTK_HOST CCTK_DEVICE void get_Press_Seed();
    //     CCTK_HOST CCTK_DEVICE void get_Z_Seed();
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE CCTK_REAL get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq);
    //     CCTK_HOST CCTK_DEVICE void WZ2Prim();

    //     /* Destructor */
    //     CCTK_HOST CCTK_DEVICE ~idealFluid();
    // };
}
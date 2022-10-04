#include <idealFluid.hxx>

/* Constructor */

using namespace AsterX;

/* Constructor */
CCTK_HOST CCTK_DEVICE idealFluid::idealFluid(CCTK_REAL gamma,
                                             CCTK_REAL (&cons)[NCONS],
                                             CCTK_REAL (&prim)[NPRIMS],
                                             CCTK_REAL (&g_lo)[4][4],
                                             CCTK_REAL (&g_up)[4][4]) {
  GammaIdealFluid = gamma;

  PrimitiveVarsSeed[RHO] = prim[RHO];
  PrimitiveVarsSeed[V1_CON] = prim[V1_CON];
  PrimitiveVarsSeed[V2_CON] = prim[V2_CON];
  PrimitiveVarsSeed[V3_CON] = prim[V3_CON];
  PrimitiveVarsSeed[EPS] = prim[EPS];

  ConservedVars[D] = cons[D];
  ConservedVars[S1_COV] = cons[S1_COV];
  ConservedVars[S2_COV] = cons[S2_COV];
  ConservedVars[S3_COV] = cons[S3_COV];
  ConservedVars[EPS] = cons[EPS];
  ConservedVars[B1] = cons[B1];
  ConservedVars[B2] = cons[B2];
  ConservedVars[B3] = cons[B3];

  gcov[TT] = g_lo[0][0];
  gcov[TX] = g_lo[0][1];
  gcov[TY] = g_lo[0][2];
  gcov[TZ] = g_lo[0][3];
  gcov[XX] = g_lo[1][1];
  gcov[XY] = g_lo[1][2];
  gcov[XZ] = g_lo[1][3];
  gcov[YY] = g_lo[2][2];
  gcov[YZ] = g_lo[2][3];
  gcov[ZZ] = g_lo[3][3];

  gcon[TT] = g_up[0][0];
  gcon[TX] = g_up[0][1];
  gcon[TY] = g_up[0][2];
  gcon[TZ] = g_up[0][3];
  gcon[XX] = g_up[1][1];
  gcon[XY] = g_up[1][2];
  gcon[XZ] = g_up[1][3];
  gcon[YY] = g_up[2][2];
  gcon[YZ] = g_up[2][3];
  gcon[ZZ] = g_up[3][3];

  CCTK_REAL alp = sqrt(-1. / gcon[TT]);

  /* B^i S_i */
  BiSi = ConservedVars[B1] * ConservedVars[S1_COV] +
         ConservedVars[B2] * ConservedVars[S2_COV] +
         ConservedVars[B3] * ConservedVars[S3_COV];

  /* Seed Lorentz factor */
  // covariant Valencia velocity:
  CCTK_REAL v1_cov = gcov[XX] * PrimitiveVarsSeed[V1_CON] +
                     gcov[XY] * PrimitiveVarsSeed[V2_CON] +
                     gcov[XZ] * PrimitiveVarsSeed[V3_CON];
  CCTK_REAL v2_cov = gcov[XY] * PrimitiveVarsSeed[V1_CON] +
                     gcov[YY] * PrimitiveVarsSeed[V2_CON] +
                     gcov[YZ] * PrimitiveVarsSeed[V3_CON];
  CCTK_REAL v3_cov = gcov[XZ] * PrimitiveVarsSeed[V1_CON] +
                     gcov[YZ] * PrimitiveVarsSeed[V2_CON] +
                     gcov[ZZ] * PrimitiveVarsSeed[V3_CON];

  CCTK_REAL vsq = v1_cov * PrimitiveVarsSeed[V1_CON] +
                  v2_cov * PrimitiveVarsSeed[V2_CON] +
                  v3_cov * PrimitiveVarsSeed[V3_CON];

  if ((vsq < 0.) && (fabs(vsq) < 1.0e-13)) {
    vsq = fabs(vsq);
  }

//   if (vsq < 0. || vsq > 1.) {
//     printf(
//         "WARNING: "
//         "vsq is either less than 0.0 or greater than 1.0, having value = %f\n",
//         vsq);
//   }

  W_Seed = 1. / sqrt(1. - vsq);

  // Bsq and bsq:
  CCTK_REAL B1_cov = gcov[XX] * ConservedVars[B1] +
                     gcov[XY] * ConservedVars[B2] +
                     gcov[XZ] * ConservedVars[B3];
  CCTK_REAL B2_cov = gcov[XY] * ConservedVars[B1] +
                     gcov[YY] * ConservedVars[B2] +
                     gcov[YZ] * ConservedVars[B3];
  CCTK_REAL B3_cov = gcov[XZ] * ConservedVars[B1] +
                     gcov[YZ] * ConservedVars[B2] +
                     gcov[ZZ] * ConservedVars[B3];

  Bsq = B1_cov * ConservedVars[B1] + B2_cov * ConservedVars[B2] +
        B3_cov * ConservedVars[B3];

  CCTK_REAL bt = W_Seed / alp *
                 (ConservedVars[B1] * v1_cov + ConservedVars[B2] * v2_cov +
                  ConservedVars[B3] * v3_cov);

  bsq = (Bsq + (alp * bt) * (alp * bt)) / (W_Seed * W_Seed);

  /* update rho seed from D and gamma */
  // rho consistent with con[D] should be better guess than rho from last
  // timestep
  PrimitiveVarsSeed[RHO] = ConservedVars[D] / W_Seed;
}

/* Called by 2dNRNoble */
CCTK_HOST CCTK_DEVICE void idealFluid::get_Ssq_Exact() {

  /* calculate S_squared */
  Ssq = ConservedVars[S1_COV] *
        (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] +
         gcon[XZ] * ConservedVars[S3_COV]);
  Ssq += ConservedVars[S2_COV] *
         (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] +
          gcon[YZ] * ConservedVars[S3_COV]);
  Ssq += ConservedVars[S3_COV] *
         (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] +
          gcon[ZZ] * ConservedVars[S3_COV]);
  if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13)) {
    Ssq = fabs(Ssq);
  }
}

CCTK_HOST CCTK_DEVICE void idealFluid::get_Press_Seed() {
  Press_Seed =
      PrimitiveVarsSeed[RHO] * PrimitiveVarsSeed[EPS] * (GammaIdealFluid - 1.0);
}

CCTK_HOST CCTK_DEVICE void idealFluid::get_Z_Seed() {
  Z_Seed = (PrimitiveVarsSeed[RHO] +
            PrimitiveVarsSeed[EPS] * PrimitiveVarsSeed[RHO] + Press_Seed) *
           W_Seed * W_Seed;
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_f0(CCTK_REAL Z,
                                                             CCTK_REAL Vsq) {
  return (Vsq * (Bsq + Z) * (Bsq + Z) -
          (BiSi * BiSi * (Bsq + 2.0 * Z)) / (Z * Z) - Ssq);
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_f1(CCTK_REAL Z,
                                                             CCTK_REAL Vsq) {
  CCTK_REAL Press = get_Press_funcZVsq(Z, Vsq);
  return ConservedVars[TAU] + ConservedVars[D] - Bsq / 2.0 * (1 + Vsq) +
         BiSi * BiSi / 2.0 / (Z * Z) - Z + Press;
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_df0dZ(CCTK_REAL Z,
                                                                CCTK_REAL Vsq) {
  return (2.0 * Vsq * (Bsq + Z) - 2.0 * BiSi * BiSi / (Z * Z) +
          2.0 * BiSi * BiSi * (Bsq + 2.0 * Z) / (Z * Z * Z));
}

CCTK_HOST CCTK_DEVICE CCTK_REAL
idealFluid::get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq) {
  return (Bsq + Z) * (Bsq + Z);
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_df1dZ(CCTK_REAL Z,
                                                                CCTK_REAL Vsq) {
  return (-BiSi * BiSi / (Z * Z * Z) - 1.0 + get_dPdZ_funcZVsq(Z, Vsq));
}

CCTK_HOST CCTK_DEVICE CCTK_REAL
idealFluid::get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq) {
  return (-Bsq / 2.0 + get_dPdVsq_funcZVsq(Z, Vsq));
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_Press_funcZVsq(CCTK_REAL Z,
                                                               CCTK_REAL Vsq) {
  return ((Z * (1.0 - Vsq) - ConservedVars[D] * sqrt(1.0 - Vsq)) *
          (GammaIdealFluid - 1.0) / (GammaIdealFluid));
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_dPdZ_funcZVsq(CCTK_REAL Z,
                                                              CCTK_REAL Vsq) {
  return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_dPdVsq_funcZVsq(CCTK_REAL Z,
                                                                CCTK_REAL Vsq) {
  return ((-Z + ConservedVars[D] / (2.0 * sqrt(1.0 - Vsq))) *
          (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

CCTK_HOST CCTK_DEVICE void idealFluid::WZ2Prim() {
  CCTK_REAL W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);

  PrimitiveVars[RHO] = ConservedVars[D] / W_Sol;

  PrimitiveVars[V1_CON] =
      (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] +
       gcon[XZ] * ConservedVars[S3_COV]) /
      (Z_Sol + Bsq);
  PrimitiveVars[V1_CON] += BiSi * ConservedVars[B1] / (Z_Sol * (Z_Sol + Bsq));

  PrimitiveVars[V2_CON] =
      (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] +
       gcon[YZ] * ConservedVars[S3_COV]) /
      (Z_Sol + Bsq);
  PrimitiveVars[V2_CON] += BiSi * ConservedVars[B2] / (Z_Sol * (Z_Sol + Bsq));

  PrimitiveVars[V3_CON] =
      (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] +
       gcon[ZZ] * ConservedVars[S3_COV]) /
      (Z_Sol + Bsq);
  PrimitiveVars[V3_CON] += BiSi * ConservedVars[B3] / (Z_Sol * (Z_Sol + Bsq));

  PrimitiveVars[EPS] =
      (Z_Sol * (1. - vsq_Sol) / PrimitiveVars[RHO] - 1.0) / GammaIdealFluid;
}

/* Destructor */
CCTK_HOST CCTK_DEVICE idealFluid::~idealFluid() {
  // How to destruct properly a vector?
}
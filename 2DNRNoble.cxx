#include "con2primFactory.hxx"

/***************************************************************************
2DNRNoble C2P
------------------------------------
2D-NR Noble scheme for c2p.
Sources: Noble+2006, Section 3.1 of Siegel+2018,
NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
****************************************************************************/
template <typename typeEoS>
CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble(
    CCTK_INT max_iter, CCTK_REAL tolf,
    typeEoS &plasma)
{ // Send con2primFactory object as reference to modify it,
  // and because we can not instantiate abstract class

    /* get Lorentz factor seed, calculated by constructor */
    CCTK_REAL W = plasma.W_Seed;

    /* get Ssq from cons (exact) */
    plasma.get_Ssq_Exact();

    /* get pressure seed */
    plasma.get_Press_Seed();

    /* get Z seed */
    plasma.get_Z_Seed();

    /* initialize unknowns for c2p, Z and vsq: */
    CCTK_REAL x[2];
    CCTK_REAL x_old[2];
    x[0] = fabs(plasma.Z_Seed);
    x[1] = (-1.0 + W * W) / (W * W);

    /* initialize old values */
    x_old[0] = x[0];
    x_old[1] = x[1];

    /* Start Recovery with 2D NR Solver */
    const CCTK_INT n = 2;
    const CCTK_REAL dv = (1. - 1.e-15);
    CCTK_REAL fvec[n];
    CCTK_REAL dx[n];
    CCTK_REAL fjac[n][n];

    CCTK_REAL detjac_inv;
    CCTK_REAL errf;
    plasma.Failed_2DNRNoble = 1;
    CCTK_INT k;
    for (k = 1; k <= max_iter; k++)
    {
        fvec[0] = plasma.get_2DNRNoble_f0(x[0], x[1]);
        fvec[1] = plasma.get_2DNRNoble_f1(x[0], x[1]);
        fjac[0][0] = plasma.get_2DNRNoble_df0dZ(x[0], x[1]);
        fjac[0][1] = plasma.get_2DNRNoble_df0dVsq(x[0], x[1]);
        fjac[1][0] = plasma.get_2DNRNoble_df1dZ(x[0], x[1]);
        fjac[1][1] = plasma.get_2DNRNoble_df1dVsq(x[0], x[1]);
        detjac_inv = 1.0 / (fjac[0][0] * fjac[1][1] - fjac[0][1] * fjac[1][0]);
        dx[0] = -detjac_inv * (fjac[1][1] * fvec[0] - fjac[0][1] * fvec[1]);
        dx[1] = -detjac_inv * (-fjac[1][0] * fvec[0] + fjac[0][0] * fvec[1]);

        errf = 0.0;
        for (CCTK_INT i = 0; i < n; i++)
        {
            errf += fabs(fvec[i]);
        }
        if (errf <= tolf)
        {
            plasma.Failed_2DNRNoble = 0;
            break;
        }

        /* save old values before calculating the new */
        x_old[0] = x[0];
        x_old[1] = x[1];

        for (CCTK_INT i = 0; i < n; i++)
        {
            x[i] += dx[i];
        }
    }
    plasma.Nit_2DNRNoble = k;

    /* Calculate primitives from Z and W */
    plasma.Z_Sol = x[0];
    plasma.vsq_Sol = x[1];
}

template CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble<AsterX::idealFluid>(
    CCTK_INT max_iter, CCTK_REAL tolf,
    AsterX::idealFluid &plasma);
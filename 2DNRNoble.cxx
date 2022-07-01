#include "con2primFactory.hxx"

/***************************************************************************
2DNRNoble C2P
------------------------------------
2D-NR Noble scheme for c2p.
Sources: Noble+2006, Section 3.1 of Siegel+2018, 
NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
****************************************************************************/
void Con2Prim_2DNRNoble(int max_iter, con2primFactory &plasma)
{ // Send con2primFactory as reference to modify it, and because we can not instantiate abstract class

    /* get Lorentz factor seed */
    plasma.get_LorentzFactor_Seed();
    double W = plasma.W_Seed;

    /* get Ssq from cons (exact) */
    plasma.get_Ssq_Exact();

    /* update rho seed from D and gamma */
    // rho consistent with con[D] should be better guess than rho from last timestep
    if (plasma.PrimitiveLabels[0] == PrimitiveLabel::P_RHO)
    {
        plasma.PrimitiveVarsSeed[0] = plasma.ConservedVars[D] / W;
    }

    /* get pressure seed */
    plasma.get_Press_Seed();

    /* get Z seed */
    plasma.get_Z_Seed();

    /* initialize unknowns for c2p, Z and vsq: */
    std::vector<double> x;
    x.reserve(2);
    x[0] = fabs(plasma.Z_Seed);
    x[1] = (-1.0 + W * W) / (W * W);

    /* Start Recovery with 2D NR Solver */
    const int n = 2;
    double tolf = 1E-15; // TODO: Make parameter
    double fvec[n];
    double dx[n];
    double fjac[n][n];

    double rho, Press, dPdx0, dPdx1;
    double detjac_inv;
    double errf;

    plasma.Failed_2DNRNoble = 1;
    for (int k = 1; k <= max_iter; k++)
    {
        Press = plasma.get_Press_funcZVsq(x[0], x[1]);
        rho = plasma.ConservedVars[D] * sqrt(1.0 - x[1]);
        fvec[0] = x[1] * x[0] * x[0] - plasma.Ssq;
        fvec[1] = plasma.ConservedVars[TAU] + plasma.ConservedVars[D] - x[0] + Press;
        dPdx0 = plasma.dPdZ_funcZVsq(x[0], x[1]);
        dPdx1 = plasma.dPdVsq_funcZVsq(x[0], x[1]);
        fjac[0][0] = 2.0 * x[1] * x[0];
        fjac[0][1] = x[0] * x[0];
        fjac[1][0] = -1.0 + dPdx0;
        fjac[1][1] = dPdx1;
        detjac_inv = 1.0 / (fjac[0][0] * fjac[1][1] - fjac[0][1] * fjac[1][0]);
        dx[0] = -detjac_inv * (fjac[1][1] * fvec[0] - fjac[0][1] * fvec[1]);
        dx[1] = -detjac_inv * (-fjac[1][0] * fvec[0] + fjac[0][0] * fvec[1]);

        errf = 0.0;
        for (int i = 0; i < n; i++)
        {
            errf += fabs(fvec[i]);
        }
        if (errf <= tolf)
        {
            plasma.Failed_2DNRNoble = 0;
            break;
        }
        for (int i = 0; i < n; i++)
        {
            x[i] += dx[i];
        }
    }

    /* Calculate primitives from Z and W */
    plasma.Z_Sol = x[0];
    plasma.vsq_Sol = x[1];
    plasma.WZ2Prim();
}

#include "con2primFactory.hxx"

/***************************************************************************
2DNRNoble C2P
------------------------------------
2D-NR Noble scheme for c2p.
Sources: Noble+2006, Section 3.1 of Siegel+2018, 
NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
****************************************************************************/
void Con2Prim_2DNRNoble(con2primFactory &plasma)
{ // Send con2primFactory object as reference to modify it, and because we can not instantiate abstract class

    /* get Lorentz factor seed */
    plasma.get_LorentzFactor_Seed();
    double W = plasma.W_Seed;

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
    double tolf = 1E-10; // TODO: Make parameter
    double fvec[n];
    double dx[n];
    double fjac[n][n];

    double Press, dPdx0, dPdx1;
    double detjac_inv;
    double errf;

    plasma.Failed_2DNRNoble = 1;
    for (int k = 1; k <= plasma.max_iterations; k++)
    {
        Press = plasma.get_Press_funcZVsq(x[0], x[1]);
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
    // Check p_m/p
    //std::cout<< plasma.bsq/Press/2.0 <<std::endl;

    /* Calculate primitives from Z and W */
    plasma.Z_Sol = x[0];
    plasma.vsq_Sol = x[1];
    plasma.WZ2Prim();
}

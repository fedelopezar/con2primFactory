#include <AsterX_Flat_idealFluid_Magnetized.hxx>

/* Constructor */

AsterX_Flat_idealFluid_Magnetized::AsterX_Flat_idealFluid_Magnetized(std::vector<double> cons, std::vector<double> prim)
{
    // Allocate some memory:
    PrimitiveLabels.reserve(8);
    PrimitiveVarsSeed.reserve(8);
    ConservedVars.reserve(8);
    PrimitiveVars.reserve(8);
    gcov.reserve(10);
    gcon.reserve(10);

    // Initialize
    get_PrimList();
    get_gcov();
    get_gcon();
    get_PrimitiveVars_Seed(prim);
    get_ConservedVars(cons, prim);
    get_Ssq_Exact();
    get_Bsq_Exact();
    get_BiSi_Exact();
}

/* Called by Constructor */

void AsterX_Flat_idealFluid_Magnetized::get_PrimList()
{
    // AsterX uses the following primitives:
    PrimitiveLabels[RHO] = PrimitiveLabel::P_RHO;
    PrimitiveLabels[V1_CON] = PrimitiveLabel::P_V1coord;
    PrimitiveLabels[V2_CON] = PrimitiveLabel::P_V2coord;
    PrimitiveLabels[V3_CON] = PrimitiveLabel::P_V3coord;
    PrimitiveLabels[EPS] = PrimitiveLabel::P_EPS;
    PrimitiveLabels[B1] = PrimitiveLabel::P_B1;
    PrimitiveLabels[B2] = PrimitiveLabel::P_B2;
    PrimitiveLabels[B3] = PrimitiveLabel::P_B3;
}

void AsterX_Flat_idealFluid_Magnetized::get_gcov()
{
    gcov[TT] = -1.0;
    gcov[TX] = 0.0;
    gcov[TY] = 0.0;
    gcov[TZ] = 0.0;
    gcov[XX] = 1.0;
    gcov[XY] = 0.0;
    gcov[XZ] = 0.0;
    gcov[YY] = 1.0;
    gcov[YZ] = 0.0;
    gcov[ZZ] = 1.0;

    sqrt_gdet = 1.0;
}

void AsterX_Flat_idealFluid_Magnetized::get_gcon()
{
    gcon[TT] = -1.0;
    gcon[TX] = 0.0;
    gcon[TY] = 0.0;
    gcon[TZ] = 0.0;
    gcon[XX] = 1.0;
    gcon[XY] = 0.0;
    gcon[XZ] = 0.0;
    gcon[YY] = 1.0;
    gcon[YZ] = 0.0;
    gcon[ZZ] = 1.0;
}

void AsterX_Flat_idealFluid_Magnetized::get_PrimitiveVars_Seed(std::vector<double> prim)
{
    // Perturb the solution primitives:
    PrimitiveVarsSeed[RHO] = 0.9 * prim[RHO]; // Use random number, iterator
    PrimitiveVarsSeed[V1_CON] = 0.9 * prim[V1_CON];
    PrimitiveVarsSeed[V2_CON] = 0.9 * prim[V2_CON];
    PrimitiveVarsSeed[V3_CON] = 0.9 * prim[V3_CON];
    PrimitiveVarsSeed[EPS] = 0.9 * prim[EPS];
    PrimitiveVarsSeed[B1] = 1.1 * prim[B1];
    PrimitiveVarsSeed[B2] = 1.1 * prim[B2];
    PrimitiveVarsSeed[B3] = 1.1 * prim[B3];
}

void AsterX_Flat_idealFluid_Magnetized::get_ConservedVars(std::vector<double> cons, std::vector<double> prim)
{
    if (cons.empty())
    {
        get_ConservedVarsFromPrimVector(prim);
    }
    else
    {
        ConservedVars[D] = cons[0];
        ConservedVars[S1_COV] = cons[1];
        ConservedVars[S2_COV] = cons[2];
        ConservedVars[S3_COV] = cons[3];
        ConservedVars[EPS] = cons[EPS];
        ConservedVars[B1] = cons[B1];
        ConservedVars[B2] = cons[B2];
        ConservedVars[B3] = cons[B3];
    }
}

void AsterX_Flat_idealFluid_Magnetized::get_ConservedVarsFromPrimVector(std::vector<double> prim)
{
    /* calculate lorentz factor: */
    // covariant coordinate velocity:
    double v1_coord_cov = gcov[XX] * prim[V1_CON] + gcov[XY] * prim[V2_CON] +
                          gcov[XZ] * prim[V3_CON];
    double v2_coord_cov = gcov[XY] * prim[V1_CON] + gcov[YY] * prim[V2_CON] +
                          gcov[YZ] * prim[V3_CON];
    double v3_coord_cov = gcov[XZ] * prim[V1_CON] + gcov[YZ] * prim[V2_CON] +
                          gcov[ZZ] * prim[V3_CON];

    // covariant Valencia velocity:
    double alp = sqrt(-1. / gcon[TT]);
    double v1_Valencia_cov = (v1_coord_cov) / alp;
    double v2_Valencia_cov = (v2_coord_cov) / alp;
    double v3_Valencia_cov = (v3_coord_cov) / alp;

    // calculate W
    double vsq = v1_Valencia_cov *
                 (gcon[XX] * v1_Valencia_cov + gcon[XY] * v2_Valencia_cov +
                  gcon[XZ] * v3_Valencia_cov);
    vsq += v2_Valencia_cov *
           (gcon[XY] * v1_Valencia_cov + gcon[YY] * v2_Valencia_cov +
            gcon[YZ] * v3_Valencia_cov);
    vsq += v3_Valencia_cov *
           (gcon[XZ] * v1_Valencia_cov + gcon[YZ] * v2_Valencia_cov +
            gcon[ZZ] * v3_Valencia_cov);
    if ((vsq < 0.) && (fabs(vsq) < 1.0e-13))
    {
        vsq = fabs(vsq);
    }
    double w_lorentz = 1. / sqrt(1. - vsq);

    double vlowx = gcov[XX] * prim[1] + gcov[XY] * prim[2] + gcov[XZ] * prim[3];
    double vlowy = gcov[XY] * prim[1] + gcov[YY] * prim[2] + gcov[XY] * prim[3];
    double vlowz = gcov[XZ] * prim[1] + gcov[YZ] * prim[2] + gcov[ZZ] * prim[3];

    // Calculate little b:
    double B1_cov = gcov[XX] * prim[B1] + gcov[XY] * prim[B2] + gcov[XZ] * prim[B3];
    double B2_cov = gcov[XY] * prim[B1] + gcov[YY] * prim[B2] + gcov[YZ] * prim[B3];
    double B3_cov = gcov[XZ] * prim[B1] + gcov[YZ] * prim[B2] + gcov[ZZ] * prim[B3];

    double bt = w_lorentz / alp * (prim[B1] * v1_Valencia_cov + prim[B2] * v2_Valencia_cov + prim[B3] * v3_Valencia_cov);
    double bx_cov = B1_cov / w_lorentz + bt * alp * v1_Valencia_cov;
    double by_cov = B2_cov / w_lorentz + bt * alp * v2_Valencia_cov;
    double bz_cov = B3_cov / w_lorentz + bt * alp * v3_Valencia_cov;

    bsq = (Bsq + (alp * bt) * (alp * bt)) / (w_lorentz * w_lorentz);

    ConservedVars[D] = sqrt_gdet * prim[RHO] * w_lorentz;

    ConservedVars[S1_COV] = sqrt_gdet * (w_lorentz * w_lorentz * (prim[RHO] + prim[RHO] * prim[EPS] * GammaIdealFluid + bsq) * v1_Valencia_cov - alp * bt * bx_cov);

    ConservedVars[S2_COV] = sqrt_gdet * (w_lorentz * w_lorentz * (prim[RHO] + prim[RHO] * prim[EPS] * GammaIdealFluid + bsq) * v2_Valencia_cov - alp * bt * by_cov);

    ConservedVars[S3_COV] = sqrt_gdet * (w_lorentz * w_lorentz * (prim[RHO] + prim[RHO] * prim[EPS] * GammaIdealFluid + bsq) * v3_Valencia_cov - alp * bt * bz_cov);

    ConservedVars[TAU] = sqrt_gdet * (w_lorentz * w_lorentz * (prim[RHO] + prim[RHO] * prim[EPS] * GammaIdealFluid + bsq) - (prim[EPS] * prim[RHO] * (GammaIdealFluid - 1.0) + bsq / 2.0) - alp * alp * bt * bt - ConservedVars[D]);

    ConservedVars[B1] = prim[B1];
    ConservedVars[B2] = prim[B2];
    ConservedVars[B3] = prim[B3];
}

/* Destructor */

AsterX_Flat_idealFluid_Magnetized::~AsterX_Flat_idealFluid_Magnetized()
{
}

// ===========
//
// Called by 2DNRNoble:
// -----------------
// ===========

void AsterX_Flat_idealFluid_Magnetized::get_LorentzFactor_Seed()
{

    // covariant coordinate velocity:
    double v1_coord_cov = gcov[XX] * PrimitiveVarsSeed[V1_CON] + gcov[XY] * PrimitiveVarsSeed[V2_CON] +
                          gcov[XZ] * PrimitiveVarsSeed[V3_CON];
    double v2_coord_cov = gcov[XY] * PrimitiveVarsSeed[V1_CON] + gcov[YY] * PrimitiveVarsSeed[V2_CON] +
                          gcov[YZ] * PrimitiveVarsSeed[V3_CON];
    double v3_coord_cov = gcov[XZ] * PrimitiveVarsSeed[V1_CON] + gcov[YZ] * PrimitiveVarsSeed[V2_CON] +
                          gcov[ZZ] * PrimitiveVarsSeed[V3_CON];

    // covariant Valencia velocity:
    double alp = sqrt(-1. / gcon[TT]);
    double v1_Valencia_cov = v1_coord_cov / alp;
    double v2_Valencia_cov = v2_coord_cov / alp;
    double v3_Valencia_cov = v3_coord_cov / alp;

    /* calculate W */
    double vsq = v1_Valencia_cov *
                 (gcon[XX] * v1_Valencia_cov + gcon[XY] * v2_Valencia_cov +
                  gcon[XZ] * v3_Valencia_cov);
    vsq += v2_Valencia_cov *
           (gcon[XY] * v1_Valencia_cov + gcon[YY] * v2_Valencia_cov +
            gcon[YZ] * v3_Valencia_cov);
    vsq += v3_Valencia_cov *
           (gcon[XZ] * v1_Valencia_cov + gcon[YZ] * v2_Valencia_cov +
            gcon[ZZ] * v3_Valencia_cov);
    if ((vsq < 0.) && (fabs(vsq) < 1.0e-13))
    {
        vsq = fabs(vsq);
    }

    /* Set the attribute Lorentz factor */
    W_Seed = 1. / sqrt(1. - vsq);
}

void AsterX_Flat_idealFluid_Magnetized::get_Ssq_Exact()
{

    /* calculate S_squared */
    Ssq =
        ConservedVars[S1_COV] * (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] +
                                 gcon[XZ] * ConservedVars[S3_COV]);
    Ssq +=
        ConservedVars[S2_COV] * (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] +
                                 gcon[YZ] * ConservedVars[S3_COV]);
    Ssq +=
        ConservedVars[S3_COV] * (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] +
                                 gcon[ZZ] * ConservedVars[S3_COV]);
    if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13))
    {
        Ssq = fabs(Ssq);
    }
}

void AsterX_Flat_idealFluid_Magnetized::get_Bsq_Exact()
{
    double B1_cov = gcov[XX] * PrimitiveVars[B1] + gcov[XY] * PrimitiveVars[B2] + gcov[XZ] * PrimitiveVars[B3];
    double B2_cov = gcov[XY] * PrimitiveVars[B1] + gcov[YY] * PrimitiveVars[B2] + gcov[YZ] * PrimitiveVars[B3];
    double B3_cov = gcov[XZ] * PrimitiveVars[B1] + gcov[YZ] * PrimitiveVars[B2] + gcov[ZZ] * PrimitiveVars[B3];
    Bsq = B1_cov * PrimitiveVars[B1] + B2_cov * PrimitiveVars[B2] + B3_cov * PrimitiveVars[B3];
}

void AsterX_Flat_idealFluid_Magnetized::get_BiSi_Exact()
{
    BiSi = PrimitiveVars[B1] * ConservedVars[S1_COV] + PrimitiveVars[B2] * ConservedVars[S2_COV] + PrimitiveVars[B3] * ConservedVars[S3_COV];
}

void AsterX_Flat_idealFluid_Magnetized::get_Press_Seed()
{
    Press_Seed = PrimitiveVarsSeed[RHO] * PrimitiveVarsSeed[EPS] * (GammaIdealFluid - 1.0);
}

void AsterX_Flat_idealFluid_Magnetized::get_Z_Seed()
{
    Z_Seed = (PrimitiveVarsSeed[RHO] + PrimitiveVarsSeed[EPS] * PrimitiveVarsSeed[RHO] + Press_Seed) * W_Seed * W_Seed;
}

double AsterX_Flat_idealFluid_Magnetized::get_2DNRNoble_f0(double Z, double Vsq)
{
    return (Vsq * (Bsq + Z) * (Bsq + Z) - (BiSi * BiSi * (Bsq + 2.0 * Z)) / (Z * Z) - Ssq);
}

double AsterX_Flat_idealFluid_Magnetized::get_2DNRNoble_f1(double Z, double Vsq)
{
    double Press = get_Press_funcZVsq(Z, Vsq);
    return ConservedVars[TAU] + ConservedVars[D] - Bsq / 2.0 * (1 + Vsq) + BiSi * BiSi / 2.0 / (Z * Z) - Z + Press;
}

double AsterX_Flat_idealFluid_Magnetized::get_2DNRNoble_df0dZ(double Z, double Vsq)
{
    return (2.0 * Vsq * (Bsq + Z) - 2.0 * BiSi * BiSi / (Z * Z) + 2.0 * BiSi * BiSi * (Bsq + 2.0 * Z) / (Z * Z * Z));
}

double AsterX_Flat_idealFluid_Magnetized::get_2DNRNoble_df0dVsq(double Z, double Vsq)
{
    return (Bsq + Z) * (Bsq + Z);
}

double AsterX_Flat_idealFluid_Magnetized::get_2DNRNoble_df1dZ(double Z, double Vsq)
{
    return (-BiSi * BiSi / (Z * Z * Z) - 1.0 + dPdZ_funcZVsq(Z, Vsq));
}

double AsterX_Flat_idealFluid_Magnetized::get_2DNRNoble_df1dVsq(double Z, double Vsq)
{
    return (-Bsq / 2.0 + dPdVsq_funcZVsq(Z, Vsq));
}

double AsterX_Flat_idealFluid_Magnetized::get_Press_funcZVsq(double Z, double Vsq)
{
    return ((Z * (1.0 - Vsq) - ConservedVars[D] * sqrt(1.0 - Vsq)) * (GammaIdealFluid - 1.0) / (GammaIdealFluid));
}

double AsterX_Flat_idealFluid_Magnetized::dPdZ_funcZVsq(double Z, double Vsq)
{
    return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

double AsterX_Flat_idealFluid_Magnetized::dPdVsq_funcZVsq(double Z, double Vsq)
{
    return ((-Z + ConservedVars[D] / (2.0 * sqrt(1.0 - Vsq))) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

void AsterX_Flat_idealFluid_Magnetized::WZ2Prim()
{
    double W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);
    // std::cout << W_Sol << std::endl;
    PrimitiveVars[RHO] = ConservedVars[D] / W_Sol;
    double alp = sqrt(-1. / gcon[TT]);

    double v1_Valencia = (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] + gcon[XZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    v1_Valencia += BiSi * PrimitiveVars[B1] / (Z_Sol * (Z_Sol + Bsq));

    double v2_Valencia = (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] + gcon[YZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    v2_Valencia += BiSi * PrimitiveVars[B2] / (Z_Sol * (Z_Sol + Bsq));

    double v3_Valencia = (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] + gcon[ZZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    v3_Valencia += BiSi * PrimitiveVars[B3] / (Z_Sol * (Z_Sol + Bsq));

    PrimitiveVars[V1_CON] = alp * v1_Valencia - alp * alp * gcon[TX];
    PrimitiveVars[V2_CON] = alp * v2_Valencia - alp * alp * gcon[TY];
    PrimitiveVars[V3_CON] = alp * v3_Valencia - alp * alp * gcon[TZ];

    PrimitiveVars[EPS] = (Z_Sol * (1. - vsq_Sol) / PrimitiveVars[RHO] - 1.0) / GammaIdealFluid;
    PrimitiveVars[B1] = ConservedVars[B1];
    PrimitiveVars[B2] = ConservedVars[B2];
    PrimitiveVars[B3] = ConservedVars[B3];
}

// ===========
//
// Called by 1DBrentPalenzuela:
// -----------------
// ===========

double AsterX_Flat_idealFluid_Magnetized::get_Press_funcRhoEps(double &rho_loc, double &eps_loc)
{
    return rho_loc * eps_loc * (GammaIdealFluid - 1.0);
}

void AsterX_Flat_idealFluid_Magnetized::xPalenzuelaToPrim()
{

    const double qPalenzuela = ConservedVars[TAU] / ConservedVars[D];
    const double rPalenzuela = Ssq / pow(ConservedVars[D], 2);
    const double sPalenzuela = Bsq / ConservedVars[D];
    const double tPalenzuela = BiSi / pow(ConservedVars[D], 3. / 2.);

    // (i)
    double Wminus2 = 1.0 - (xPalenzuela_Sol * xPalenzuela_Sol * rPalenzuela + (2 * xPalenzuela_Sol + sPalenzuela) * tPalenzuela * tPalenzuela) / (xPalenzuela_Sol * xPalenzuela_Sol * (xPalenzuela_Sol + sPalenzuela) * (xPalenzuela_Sol + sPalenzuela));
    Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
    const double W_loc = pow(Wminus2, -0.5);

    // (ii)
    PrimitiveVars[RHO] = ConservedVars[D] / W_loc;

    // (iii)
    PrimitiveVars[EPS] = W_loc - 1.0 + (1.0 - W_loc * W_loc) * xPalenzuela_Sol / W_loc + W_loc * (qPalenzuela - sPalenzuela + tPalenzuela * tPalenzuela / (2 * xPalenzuela_Sol * xPalenzuela_Sol) + sPalenzuela / (2 * W_loc * W_loc));

    // (iv)
    // double P_loc = get_Press_funcRhoEps(rho_loc, eps_loc);

    // Taken from WZ2Prim (2DNRNoble)
    Z_Sol = xPalenzuela_Sol * PrimitiveVars[RHO] * W_loc;
    double alp = sqrt(-1. / gcon[TT]);

    double v1_Valencia = (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] + gcon[XZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    v1_Valencia += BiSi * PrimitiveVars[B1] / (Z_Sol * (Z_Sol + Bsq));

    double v2_Valencia = (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] + gcon[YZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    v2_Valencia += BiSi * PrimitiveVars[B2] / (Z_Sol * (Z_Sol + Bsq));

    double v3_Valencia = (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] + gcon[ZZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    v3_Valencia += BiSi * PrimitiveVars[B3] / (Z_Sol * (Z_Sol + Bsq));

    PrimitiveVars[V1_CON] = alp * v1_Valencia - alp * alp * gcon[TX];
    PrimitiveVars[V2_CON] = alp * v2_Valencia - alp * alp * gcon[TY];
    PrimitiveVars[V3_CON] = alp * v3_Valencia - alp * alp * gcon[TZ];

    PrimitiveVars[B1] = ConservedVars[B1];
    PrimitiveVars[B2] = ConservedVars[B2];
    PrimitiveVars[B3] = ConservedVars[B3];
}
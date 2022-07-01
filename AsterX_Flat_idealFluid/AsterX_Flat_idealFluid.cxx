#include <AsterX_Flat_idealFluid.hxx>

/* Constructor */

AsterX_Flat_idealFluid::AsterX_Flat_idealFluid(std::vector<double> cons, std::vector<double> prim)
{
    // Allocate some memory:
    PrimitiveLabels.reserve(5);
    PrimitiveVarsSeed.reserve(5);
    ConservedVars.reserve(5);
    PrimitiveVars.reserve(5);
    gcov.reserve(10);
    gcon.reserve(10);

    // Initialize
    get_PrimList();
    get_gcov();
    get_gcon();
    get_PrimitiveVars_Seed(prim);
    get_ConservedVars(cons, prim);
}

/* Called by Constructor */

void AsterX_Flat_idealFluid::get_PrimList()
{
    // AsterX uses the following primitives:
    PrimitiveLabels[RHO] = PrimitiveLabel::P_RHO;
    PrimitiveLabels[V1_CON] = PrimitiveLabel::P_EPS;
    PrimitiveLabels[V2_CON] = PrimitiveLabel::P_V1coord;
    PrimitiveLabels[V3_CON] = PrimitiveLabel::P_V2coord;
    PrimitiveLabels[EPS] = PrimitiveLabel::P_V3coord;
}

void AsterX_Flat_idealFluid::get_gcov()
{
    // This con2primFactory is restricted to flat spacetime
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

void AsterX_Flat_idealFluid::get_gcon()
{
    // This con2primFactory is restricted to flat spacetime
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

void AsterX_Flat_idealFluid::get_PrimitiveVars_Seed(std::vector<double> prim)
{
    // Perturb the solution primitives:
    PrimitiveVarsSeed[RHO] = 0.9 * prim[RHO]; // Use random number, iterator
    PrimitiveVarsSeed[V1_CON] = prim[V1_CON];
    PrimitiveVarsSeed[V2_CON] = prim[V2_CON];
    PrimitiveVarsSeed[V3_CON] = prim[V3_CON];
    PrimitiveVarsSeed[EPS] = 0.9 * prim[EPS];
}

void AsterX_Flat_idealFluid::get_ConservedVars(std::vector<double> cons, std::vector<double> prim)
{
    if (cons.empty())
    {
        get_ConservedVarsFromPrimVector(prim);
    }
    else
    {
        ConservedVars[D] = cons[0]; // Use iterator
        ConservedVars[S1_COV] = cons[1];
        ConservedVars[S2_COV] = cons[2];
        ConservedVars[S3_COV] = cons[3];
        ConservedVars[4] = cons[4];
    }
}

void AsterX_Flat_idealFluid::get_ConservedVarsFromPrimVector(std::vector<double> prim)
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
    double v1_Valencia_cov = (v1_coord_cov + gcov[TX]) / alp;
    double v2_Valencia_cov = (v2_coord_cov + gcov[TY]) / alp;
    double v3_Valencia_cov = (v3_coord_cov + gcov[TZ]) / alp;

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

    ConservedVars[D] = sqrt_gdet * prim[RHO] * w_lorentz;
    ConservedVars[S1_COV] = sqrt_gdet * prim[RHO] * w_lorentz * w_lorentz * (1.0 + prim[EPS] * GammaIdealFluid) * vlowx;
    ConservedVars[S2_COV] = sqrt_gdet * prim[RHO] * w_lorentz * w_lorentz * (1.0 + prim[EPS] * GammaIdealFluid) * vlowy;
    ConservedVars[S3_COV] = sqrt_gdet * prim[RHO] * w_lorentz * w_lorentz * (1.0 + prim[EPS] * GammaIdealFluid) * vlowz;
    ConservedVars[TAU] = sqrt_gdet * prim[RHO] * w_lorentz * w_lorentz * (1.0 + prim[EPS] * GammaIdealFluid) - sqrt_gdet * prim[EPS] * prim[RHO] * (GammaIdealFluid - 1.0) - sqrt_gdet * ConservedVars[D];
}

/* Called by 2dNRNoble */

void AsterX_Flat_idealFluid::get_LorentzFactor_Seed() // From prims
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
    double v1_Valencia_cov = (v1_coord_cov + gcov[TX]) / alp;
    double v2_Valencia_cov = (v2_coord_cov + gcov[TY]) / alp;
    double v3_Valencia_cov = (v3_coord_cov + gcov[TZ]) / alp;

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

void AsterX_Flat_idealFluid::get_Ssq_Exact()
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

void AsterX_Flat_idealFluid::get_Press_Seed()
{
    Press_Seed = PrimitiveVarsSeed[RHO] * PrimitiveVarsSeed[EPS] * (GammaIdealFluid - 1.0);
}

void AsterX_Flat_idealFluid::get_Z_Seed()
{
    Z_Seed = (PrimitiveVarsSeed[RHO] + PrimitiveVarsSeed[EPS] * PrimitiveVarsSeed[RHO] + Press_Seed) * W_Seed * W_Seed;
}

double AsterX_Flat_idealFluid::get_Press_funcZVsq(double Z, double Vsq)
{
    return ((Z * (1.0 - Vsq) - ConservedVars[D] * sqrt(1.0 - Vsq)) * (GammaIdealFluid - 1.0) / (GammaIdealFluid));
};

double AsterX_Flat_idealFluid::dPdZ_funcZVsq(double Z, double Vsq)
{
    return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

double AsterX_Flat_idealFluid::dPdVsq_funcZVsq(double Z, double Vsq)
{
    return ((-Z + ConservedVars[D] / (2.0 * sqrt(1.0 - Vsq))) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

void AsterX_Flat_idealFluid::WZ2Prim()
{
    double W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);
    PrimitiveVars[RHO] = ConservedVars[D] / W_Sol;
    double v1_Valencia_cov = ConservedVars[S1_COV] / Z_Sol;
    double v2_Valencia_cov = ConservedVars[S2_COV] / Z_Sol;
    double v3_Valencia_cov = ConservedVars[S3_COV] / Z_Sol;
    double alp = sqrt(-1. / gcon[TT]);
    double v1_coord_cov = alp * v1_Valencia_cov - gcov[TX];
    double v2_coord_cov = alp * v2_Valencia_cov - gcov[TY];
    double v3_coord_cov = alp * v3_Valencia_cov - gcov[TZ];
    PrimitiveVars[V1_CON] = gcon[XX] * v1_coord_cov + gcon[XY] * v2_coord_cov +
                            gcon[XZ] * v3_coord_cov;
    PrimitiveVars[V2_CON] = gcon[XY] * v1_coord_cov + gcon[YY] * v2_coord_cov +
                            gcon[YZ] * v3_coord_cov;
    PrimitiveVars[V3_CON] = gcon[XZ] * v1_coord_cov + gcon[YZ] * v2_coord_cov +
                            gcon[ZZ] * v3_coord_cov;
    PrimitiveVars[EPS] = (Z_Sol * (1. - vsq_Sol) / PrimitiveVars[RHO] - 1.0) / GammaIdealFluid;
}

/* Destructor */

AsterX_Flat_idealFluid::~AsterX_Flat_idealFluid()
{
    // How to destruct properly a vector?
}
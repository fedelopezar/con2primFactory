#include <iostream>
#include <vector>
#include <math.h>

/* Macros for indices of vector ConservedVars */
#define D 0
#define S1_COV 1
#define S2_COV 2
#define S3_COV 3
#define TAU 4
#define B1 5
#define B2 6
#define B3 7

/* Type Primitive to define your specific primitives */
enum class PrimitiveLabel
{
    P_RHO,
    P_UU,
    P_EPS,
    P_Temp,
    P_V1Valencia,
    P_V1coord,
    P_V2Valencia,
    P_V2coord,
    P_V3Valencia,
    P_V3coord,
    P_B1,
    P_B2,
    P_B3
};

/* Abstract class con2primFactory */
class con2primFactory
{
public:
    /* The constructor must initialize the following vectors */
    std::vector<PrimitiveLabel> PrimitiveLabels;
    std::vector<double> ConservedVars;     // Conserved to solve
    std::vector<double> PrimitiveVarsSeed; // Primitive seeds
    std::vector<double> PrimitiveVars;     // Primitive solution
    std::vector<double> gcov;
    std::vector<double> gcon;
    double sqrt_gdet;
    double Bsq, bsq, BiSi, Ssq;
    int max_iterations = 1000;

    virtual void get_Ssq_Exact() = 0; // From cons (exact)
    virtual void get_Bsq_Exact() = 0; // From cons (exact)
    virtual void get_BiSi_Exact() = 0; // From cons (exact)

    /* These must be set for 2DNRNoble scheme */
    int Failed_2DNRNoble;
    double W_Seed, vsq_Sol, Press_Seed, Z_Seed, Z_Sol, vsq_Seed;
    virtual void get_LorentzFactor_Seed() = 0; // From seed prims and cons
    virtual void get_Press_Seed() = 0;         // From seed prims and cons
    virtual void get_Z_Seed() = 0;             // From seed prims and cons
    virtual double get_2DNRNoble_f0(double Z, double Vsq) = 0;
    virtual double get_2DNRNoble_f1(double Z, double Vsq) = 0;
    virtual double get_2DNRNoble_df0dZ(double Z, double Vsq) = 0;
    virtual double get_2DNRNoble_df0dVsq(double Z, double Vsq) = 0;
    virtual double get_2DNRNoble_df1dZ(double Z, double Vsq) = 0;
    virtual double get_2DNRNoble_df1dVsq(double Z, double Vsq) = 0;

    virtual double get_Press_funcZVsq(double Z, double Vsq) = 0;
    virtual double dPdZ_funcZVsq(double Z, double Vsq) = 0;
    virtual double dPdVsq_funcZVsq(double Z, double Vsq) = 0;
    virtual void WZ2Prim() = 0;

    /* Called by 1DBrentPalenzuela */
    int Failed_1DBrentPalenzuela;
    double xPalenzuela_Sol;
    virtual double get_Press_funcRhoEps(double &rho_loc, double &eps_loc) = 0;
    virtual void xPalenzuelaToPrim() = 0;
};

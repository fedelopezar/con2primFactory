#include "con2primFactory.hxx"

/* Macros for indices of vector PrimitiveVars */
#define RHO 0
#define V1_COV 1
#define V1_CON 1
#define V2_COV 2
#define V2_CON 2
#define V3_COV 3
#define V3_CON 3
#define EPS 4

/* Macros for indices of gcov/gcon */
#define TT 0
#define TX 1
#define TY 2
#define TZ 3
#define XX 4
#define XY 5
#define XZ 6
#define YY 7
#define YZ 8
#define ZZ 9

class AsterX_Flat_idealFluid : public con2primFactory
{
public:
    /* Some attributes */
    double GammaIdealFluid = 4. / 3.;
    std::vector<double> gcov;
    std::vector<double> gcon;
    double sqrt_gdet;

    /* Constructor */
    AsterX_Flat_idealFluid(std::vector<double> cons, std::vector<double> prim);
    /* Called by constructor */
    void get_PrimList();
    void get_gcov();
    void get_gcon();
    void get_PrimitiveVars_Seed(std::vector<double> prim);
    void get_ConservedVars(std::vector<double> cons, std::vector<double> prim);
    void get_ConservedVarsFromPrimVector(std::vector<double> prim);

    /* Called by 2DNRNoble */
    void get_LorentzFactor_Seed();
    void get_Ssq_Exact();
    void get_Press_Seed();
    void get_Z_Seed();
    double get_Press_funcZVsq(double Z, double Vsq);
    double dPdZ_funcZVsq(double Z, double Vsq);
    double dPdVsq_funcZVsq(double Z, double Vsq);
    void WZ2Prim();

    /* Destructor */
    ~AsterX_Flat_idealFluid();
};

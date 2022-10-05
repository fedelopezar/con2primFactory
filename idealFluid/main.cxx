#include "idealFluid.hxx"

#define CCTK_ARGUMENTS
#define CCTK_PASS_CTOC

// For testing purposes:
#include "H5Cpp.h"
const H5std_string FILE_NAME("2DNRNoble.h5");
const H5std_string DATASET0_NAME("prims_or");
const H5std_string DATASET1_NAME("prims_sol");
const H5std_string DATASET2_NAME("N_RHO");
const H5std_string DATASET3_NAME("N_TEMP");
//

using namespace AsterX;

extern template void Con2Prim_2DNRNoble<idealFluid>(
    CCTK_INT max_iter, CCTK_REAL tolf,
    idealFluid &plasma);

template <typename typeEoS>
int Test_Con2Prim_typeEoS(CCTK_ARGUMENTS)
{

    CCTK_INT max_iter = 100;
    CCTK_REAL c2p_tol = 1E-10;
    CCTK_REAL gamma = 4. / 3.;

    CCTK_REAL conserved_0[NCONS] = {0};
    CCTK_REAL prims_0[NPRIMS] = {0};
    prims_0[1] = 0.8;
    //prims_0[7] = 0.001;

    CCTK_REAL g_up[4][4] = {0};
    CCTK_REAL g_lo[4][4] = {0};
    g_lo[0][0] = -1.0;
    g_lo[1][1] = 1.0;
    g_lo[2][2] = 1.0;
    g_lo[3][3] = 1.0;
    g_up[0][0] = -1.0;
    g_up[1][1] = 1.0;
    g_up[2][2] = 1.0;
    g_up[3][3] = 1.0;

    CCTK_REAL cons[NCONS] = {0};
    CCTK_REAL prims[NPRIMS] = {0};

    /* Parameters for the test */
    double RHO_CGS2CCTK = 1.6189988336901327e-18;
    double RHO_MIN = 1.0E6 * RHO_CGS2CCTK;
    double RHO_MAX = 1.0E14 * RHO_CGS2CCTK;
    int N_RHO = 10;
    double RHO_EXP_STEP = (log10(RHO_MAX) - log10(RHO_MIN)) / N_RHO;
    double RHO_EXP = log10(RHO_MIN);

    double TEMP_MIN = 1.0E8;
    double TEMP_MAX = 1.0E10;
    int N_TEMP = 15;
    double TEMP_EXP_STEP = (log10(TEMP_MAX) - log10(TEMP_MIN)) / N_TEMP;
    double TEMP_EXP = log10(TEMP_MIN);
    double PRESS_CGS2CCTK = 1.8013791430560688e-39;

    /* Loop over RHO and TEMP */
    double rho;
    double temp, eps;

    double prims_or[N_RHO * N_TEMP][8];
    double prims_sol[N_RHO * N_TEMP][8];
    int count = 0;
    double B;
    for (int i = 1; i <= N_RHO; i++)
    {
        // rho = pow(10.0, RHO_EXP) * RHO_CGS2CCTK;
        rho = pow(10.0, RHO_EXP);
        RHO_EXP += RHO_EXP_STEP;

        TEMP_EXP = log10(TEMP_MIN);
        for (int j = 1; j <= N_TEMP; j++)
        {
            temp = pow(10.0, TEMP_EXP);
            TEMP_EXP += TEMP_EXP_STEP;

            eps = temp / (gamma - 1.0);
            prims_0[0] = rho;
            prims_0[4] = eps;
            prims_0[7] = sqrt(eps * rho) / 10;
            prims_or[count][0] = rho;
            prims_or[count][1] = prims_0[1];
            prims_or[count][2] = prims_0[2];
            prims_or[count][3] = prims_0[3];
            prims_or[count][4] = eps;
            prims_or[count][5] = prims_0[5];
            prims_or[count][6] = prims_0[6];
            prims_or[count][7] = prims_0[7];

            // Construct con2primFactory object:
            typeEoS plasma_0(gamma, conserved_0, prims_0, g_lo, g_up, 1);

            // 1) Try 2DNRNoble
            Con2Prim_2DNRNoble(max_iter, c2p_tol, plasma_0);

            plasma_0.WZ2Prim();

            prims_sol[count][0] = plasma_0.PrimitiveVars[0];
            prims_sol[count][1] = plasma_0.PrimitiveVars[1];
            prims_sol[count][2] = plasma_0.PrimitiveVars[2];
            prims_sol[count][3] = plasma_0.PrimitiveVars[3];
            prims_sol[count][4] = plasma_0.PrimitiveVars[4];
            prims_sol[count][5] = prims_0[5];
            prims_sol[count][6] = prims_0[6];
            prims_sol[count][7] = prims_0[7];
            count++;
        }
    }

    /* Write test results */
    {
        // Try block to detect exceptions raised by any of the calls inside it
        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            H5::Exception::dontPrint();

            // Create a new file using the default property lists.
            H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);

            // Create the data space for the dataset.
            hsize_t dims[2]; // dataset dimensions
            dims[0] = N_RHO * N_TEMP;
            dims[1] = 8;
            hsize_t dims_scalar[1]; // dataset dimensions
            dims_scalar[0] = 1;
            int dataset_scalar[1];
            const int RANK = 2;
            H5::DataSpace dataspace0(RANK, dims); // TODO: Create group
            H5::DataSpace dataspace1(RANK, dims);
            H5::DataSpace dataspace2(1, dims_scalar);
            H5::DataSpace dataspace3(1, dims_scalar);

            // Create the dataset.
            H5::DataSet dataset0 = file.createDataSet(DATASET0_NAME, H5::PredType::NATIVE_DOUBLE, dataspace0);
            H5::DataSet dataset1 = file.createDataSet(DATASET1_NAME, H5::PredType::NATIVE_DOUBLE, dataspace1);
            H5::DataSet dataset2 = file.createDataSet(DATASET2_NAME, H5::PredType::NATIVE_INT, dataspace2);
            H5::DataSet dataset3 = file.createDataSet(DATASET3_NAME, H5::PredType::NATIVE_INT, dataspace3);

            dataset0.write(prims_or, H5::PredType::NATIVE_DOUBLE);
            dataset1.write(prims_sol, H5::PredType::NATIVE_DOUBLE);
            dataset_scalar[0] = N_RHO;
            dataset2.write(dataset_scalar, H5::PredType::NATIVE_INT);
            dataset_scalar[0] = N_TEMP;
            dataset3.write(dataset_scalar, H5::PredType::NATIVE_INT);

        } // end of try block

        // catch failure caused by the H5File operations
        catch (H5::FileIException error)
        {
            error.printErrorStack();
            return -1;
        }

        // catch failure caused by the DataSet operations
        catch (H5::DataSetIException error)
        {
            error.printErrorStack();
            return -1;
        }

        // catch failure caused by the DataSpace operations
        catch (H5::DataSpaceIException error)
        {
            error.printErrorStack();
            return -1;
        }
    }

    return 0;
}

/***************************************************************************
 * AsterX_Con2Prim
 * ---
 *  Routines implemented:
 *   1) 2DNRNoble
 *
 *   Based on con2primFactory (https://github.com/fedelopezar/con2primFactory)
 *   ****************************************************************************/

template <typename typeEoS>
void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS)
{
    CCTK_INT max_iter = 100;
    CCTK_REAL c2p_tol = 1E-10;
    CCTK_REAL gamma = 4. / 3.;

    CCTK_REAL g_up[4][4];
    CCTK_REAL g_lo[4][4];
    CCTK_REAL cons[NCONS];

    cons[D] = 0.0;
    cons[S1_COV] = 0.0;
    cons[S2_COV] = 0.0;
    cons[S3_COV] = 0.0;
    cons[TAU] = 0.0;
    cons[B1] = 0.0;
    cons[B2] = 0.0;
    cons[B3] = 0.0;

    CCTK_REAL prims[NPRIMS];
    prims[RHO] = 0.0;
    prims[V1_CON] = 0.0;
    prims[V2_CON] = 0.0;
    prims[V3_CON] = 0.0;
    prims[EPS] = 0.0;

    // Construct con2primFactory object:
    typeEoS plasma_0(gamma, cons, prims, g_lo, g_up, 0);
    // 1) Try 2DNRNoble
    Con2Prim_2DNRNoble(max_iter, c2p_tol, plasma_0);
}

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS)
{
    if (1)
    { // Use this if for idealFluid/tabeos
        AsterX_Con2Prim_typeEoS<idealFluid>(CCTK_PASS_CTOC);
        // CCTK_PASS_CTOC == cctkGH, and more. Preferred over just cctkGH.
    }
}

int main(int argn, char **argv)
{
    int testFailed = Test_Con2Prim_typeEoS<idealFluid>();
    return 0;
}

#include "AsterX_Flat_idealFluid_Magnetized.hxx"
#include "H5Cpp.h"

const H5std_string FILE_NAME("2DNRNoble.h5");
const H5std_string DATASET0_NAME("prims_or");
const H5std_string DATASET1_NAME("prims_sol");
const H5std_string DATASET2_NAME("N_RHO");
const H5std_string DATASET3_NAME("N_TEMP");

extern void Con2Prim_2DNRNoble(con2primFactory &plasma);
extern void Con2Prim_1DBrentPalenzuela(con2primFactory &plasma);

int main(int argn, char **argv)
{

    std::vector<double> conserved_0;
    std::vector<double> prims_0{0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001};

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

            eps = temp / (4.0 / 3.0 - 1.0);
            prims_0[0] = rho*0.9;
            prims_0[4] = eps*0.9;
            prims_0[7] = sqrt(eps * rho) / 10;
            prims_or[count][0] = rho;
            prims_or[count][1] = prims_0[1];
            prims_or[count][2] = prims_0[2];
            prims_or[count][3] = prims_0[3];
            prims_or[count][4] = eps;
            prims_or[count][5] = prims_0[5];
            prims_or[count][6] = prims_0[6];
            prims_or[count][7] = prims_0[7];

            AsterX_Flat_idealFluid_Magnetized plasma_0(conserved_0, prims_0);
            Con2Prim_2DNRNoble(plasma_0);
            Con2Prim_1DBrentPalenzuela(plasma_0);

            std::cout << "0: " << plasma_0.PrimitiveVars[0] << " " << prims_or[count][0] << "\n";
            std::cout << "1: " << plasma_0.PrimitiveVars[1] << " " << prims_or[count][1] << "\n";
            std::cout << "2: " << plasma_0.PrimitiveVars[2] << " " << prims_or[count][2] << "\n";
            std::cout << "3: " << plasma_0.PrimitiveVars[3] << " " << prims_or[count][3] << "\n";
            std::cout << "4: " << plasma_0.PrimitiveVars[4] << " " << prims_or[count][4] << "\n";

            prims_sol[count][0] = plasma_0.PrimitiveVars[0];
            prims_sol[count][1] = plasma_0.PrimitiveVars[1];
            prims_sol[count][2] = plasma_0.PrimitiveVars[2];
            prims_sol[count][3] = plasma_0.PrimitiveVars[3];
            prims_sol[count][4] = plasma_0.PrimitiveVars[4];
            prims_sol[count][5] = plasma_0.PrimitiveVars[5];
            prims_sol[count][6] = plasma_0.PrimitiveVars[6];
            prims_sol[count][7] = plasma_0.PrimitiveVars[7];
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
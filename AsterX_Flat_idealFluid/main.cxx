#include "AsterX_Flat_idealFluid.hxx"
extern void Con2Prim_2DNRNoble(int max_iter, con2primFactory& plasma);
int main(int argn, char **argv)
{
    std::vector<double> conserved_0;
    std::vector<double> prims_0{1.0, .0, 0.0, 0.0, 0.1};
    std::cout << prims_0[0] << std::endl;
    std::cout << prims_0[1] << std::endl;
    std::cout << prims_0[2] << std::endl;
    std::cout << prims_0[3] << std::endl;
    std::cout << prims_0[4] << std::endl;
    AsterX_Flat_idealFluid plasma_0(conserved_0, prims_0);
    Con2Prim_2DNRNoble( 1000, plasma_0);


    std::cout << plasma_0.PrimitiveVars[0] << std::endl;
    std::cout << plasma_0.PrimitiveVars[1] << std::endl;
    std::cout << plasma_0.PrimitiveVars[2] << std::endl;
    std::cout << plasma_0.PrimitiveVars[3] << std::endl;
    std::cout << plasma_0.PrimitiveVars[4] << std::endl;
    return 0;
}
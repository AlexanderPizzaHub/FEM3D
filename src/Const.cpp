#include "Const.hpp"

namespace constants
{
    inline namespace test1
    {
        PetscScalar Source(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return 3*(PI*PI)*sin(PI*x)*sin(PI*y)*sin(PI*z);
        };

        PetscScalar BdryDirichlet(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return sin(PI*x)*sin(PI*y)*sin(PI*z);
        };

        PetscScalar Exact(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return sin(PI*x)*sin(PI*y)*sin(PI*z);
        }
    }
}
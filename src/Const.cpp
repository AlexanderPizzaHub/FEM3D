#include "Const.hpp"

namespace constants
{
    namespace test1
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

    inline namespace test2
    {
        PetscScalar Source(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return (y*y*y)/6;
        };

        PetscScalar BdryDirichlet(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return y;
        };

        PetscScalar Exact(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return y;
        }

        PetscScalar BdryNeumann(PetscScalar x, PetscScalar y, PetscScalar z)
        {
            return 0;
        }
    }
}
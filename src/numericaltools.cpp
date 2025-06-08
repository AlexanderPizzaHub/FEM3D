#include "numericaltools.hpp"

namespace numerical
{
    PetscErrorCode near(PetscScalar a, PetscScalar b)
    {
        return PetscAbs(a - b) < 1e-10;
    };
}
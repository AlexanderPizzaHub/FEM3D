#include "numericaltools.hpp"

namespace numerical
{
    PetscErrorCode VecMatVecInner(const Vec v1, const Mat M, const Vec v2, PetscScalar &result)
    {
        Vec tmp;
        PetscCall(VecDuplicate(v1, &tmp));
        PetscCall(MatMult(M, v2, tmp));
        PetscCall(VecDot(v1, tmp, &result));
        PetscCall(VecDestroy(&tmp));
        return 0;
    }

    PetscErrorCode VecErrL2(const Vec vec1, const Vec vec2, PetscScalar &err)
    {
        Vec residual;
        PetscCall(VecDuplicate(vec1, &residual));
        PetscCall(VecWAXPY(residual, -1.0, vec1, vec2));
        PetscCall(VecNorm(residual, NORM_2, &err));
        return 0;
    }
    PetscErrorCode VecErrL2Weight(const Vec vec1, const Mat M, const Vec vec2, PetscScalar &err)
    {
        Vec residual;
        PetscCall(VecDuplicate(vec1, &residual));
        PetscCall(VecWAXPY(residual, -1.0, vec1, vec2));
        PetscCall(VecMatVecInner(residual, M, residual, err));
        err = PetscSqrtReal(err);
        return 0;
    }

    PetscErrorCode VecErrL2Rel(const Vec vec1, const Vec vec2, PetscScalar &err)
    {
        Vec residual;
        PetscScalar vecnorm;

        PetscCall(VecDuplicate(vec1, &residual));

        PetscCall(VecWAXPY(residual, -1.0, vec1, vec2));

        PetscCall(VecNorm(residual, NORM_2, &err));
        PetscCall(VecNorm(vec2, NORM_2, &vecnorm));
        err /= vecnorm;
        return 0 ;
    }

    PetscErrorCode VecErrL2RelWeight(const Vec vec1, const Mat M, const Vec vec2, PetscScalar &err)
    {
        Vec residual;
        PetscScalar vecnorm;

        PetscCall(VecDuplicate(vec1, &residual));

        PetscCall(VecWAXPY(residual, -1.0, vec1, vec2));

        PetscCall(VecMatVecInner(residual, M, residual, err));
        PetscCall(VecMatVecInner(vec2, M, vec2, vecnorm));

        err = PetscSqrtReal(err);
        vecnorm = PetscSqrtReal(vecnorm);

        err /= vecnorm;
        return 0;
    }
}
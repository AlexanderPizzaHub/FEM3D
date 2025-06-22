/*
此模块负责例如矩阵向量乘，等功能编写。
*/
#pragma once 

#include <petsc.h>
#include <vector>

namespace numerical
{
    inline bool near(PetscScalar a, PetscScalar b){
           return PetscAbs(a - b) < 1e-10;
       };

    inline void determinant33(const PetscScalar* matrix, PetscScalar &det){
        // 计算3x3矩阵的行列式，矩阵matrix为row-major
        det = matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7]) -
              matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6]) +
              matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
    };

    inline void vecinner(const PetscScalar* vec1, const PetscScalar* vec2, PetscInt size, PetscScalar &result){
        // 计算两个向量的内积
        result = 0.0;
        for (PetscInt i = 0; i < size; ++i)
        {
            result += vec1[i] * vec2[i];
        }
    };

    PetscErrorCode VecMatVecInner(const Vec v1, const Mat M, const Vec v2, PetscScalar &result);

    PetscErrorCode VecErrL2(const Vec vec1, const Vec vec2, PetscScalar& err);
    PetscErrorCode VecErrL2Weight(const Vec vec1, const Mat M, const Vec vec2, PetscScalar& err);

    PetscErrorCode VecErrL2Rel(const Vec vec1, const Vec vec2, PetscScalar& err);
    PetscErrorCode VecErrL2RelWeight(const Vec vec1, const Mat M, const Vec vec2, PetscScalar& err);

    PetscErrorCode VecErrLinf(const Vec vec1, const Vec vec2, PetscScalar& err);

}
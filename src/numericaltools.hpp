/*
此模块负责例如矩阵向量乘，等功能编写。
*/

#include <petsc.h>
#include <vector>

namespace numerical
{
    PetscErrorCode near(PetscScalar a, PetscScalar b);
}
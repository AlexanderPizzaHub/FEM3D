/*
此模块负责一些便利工具
*/
#pragma once 

#include <petsc.h>

static int dim = 3;

namespace utils
{
    inline PetscErrorCode VecSetup(PetscInt n, Vec &v)
    {
        PetscCall(VecCreate(PETSC_COMM_WORLD, &v));
        PetscCall(VecSetFromOptions(v));
        PetscCall(VecSetSizes(v, n, PETSC_DECIDE));
        return 0;
    }

    inline PetscErrorCode MatSetup(PetscInt n, PetscInt m, Mat &M)
    {
        PetscCall(MatCreate(PETSC_COMM_WORLD, &M));
        PetscCall(MatSetFromOptions(M));
        PetscCall(MatSetSizes(M, n, m, PETSC_DECIDE, PETSC_DECIDE));
        return 0;
    }
}
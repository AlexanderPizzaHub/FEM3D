/*
此模块负责一些常数的编写。
*/
#pragma once 

#include <petsc.h>

namespace constants
{
    const PetscScalar PI = 3.1415926535897932384626;
    namespace test1
    {
        PetscScalar Source(PetscScalar x, PetscScalar y, PetscScalar z);
        
        PetscScalar BdryDirichlet(PetscScalar x, PetscScalar y, PetscScalar z);
        
        PetscScalar Exact(PetscScalar x, PetscScalar y, PetscScalar z);
        
    };

    inline namespace test2
    {
        PetscScalar Source(PetscScalar x, PetscScalar y, PetscScalar z);

        PetscScalar BdryDirichlet(PetscScalar x, PetscScalar y, PetscScalar z);

        PetscScalar Exact(PetscScalar x, PetscScalar y, PetscScalar z);

        PetscScalar BdryNeumann(PetscScalar x, PetscScalar y, PetscScalar z);
    }
}
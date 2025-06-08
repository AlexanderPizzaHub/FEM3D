/*
此模块负责FEM边界的处理。暂时不写。

负责给网格点打标签。
*/

#pragma once 
#include <petsc.h>
#include "FEMmachine.hpp"


namespace fempatch
{
    class FEMPatchBase
    {
        public:
            FEMPatchBase(femm::LagrangeFEMBase &fem);
            virtual ~FEMPatchBase() = default;

            virtual PetscErrorCode ApplyBC(Vec bc) = 0; // 应用边界条件

        protected:
            femm::LagrangeFEMBase *fem_;

            virtual PetscErrorCode MakePatchNodeIndices() = 0;

    };
}


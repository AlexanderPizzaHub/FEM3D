/*
此模块负责FEM边界的处理。暂时不写。

负责给网格点打标签。
*/

#pragma once 
#include <petsc.h>
#include <iostream>
#include "FEMmachine.hpp"
#include "mesh.hpp"
#include "numericaltools.hpp"
#include "utils.hpp"
#include "FEMpatch.hpp"

namespace fempatch
{
    class FEMPatchBase
    {
        public:
            FEMPatchBase(femm::LagrangeFEMBase &fem);
            virtual ~FEMPatchBase() = default;

            virtual PetscErrorCode ApplyBC(const Mat stiff_glb, Mat& stiff_inner, const Vec bc_glb, Vec& bc_inner) = 0; // 应用边界条件
            virtual PetscErrorCode AugVec(const Vec vec_original, Vec& vec_aug) = 0;

        protected:
            femm::LagrangeFEMBase *fem_;

            virtual PetscErrorCode MakePatchNodeIndices() = 0;
            

            IS inner_nodes_is_; // 内节点的索引集 // 这个不应该有
            IS patch_nodes_is_; // 补丁节点的索引集

    };

    class FEMPatchDirichletZero : public FEMPatchBase
    {
        /*
        Dirichlet边界条件，零值
        */
        public:
            FEMPatchDirichletZero(femm::LagrangeFEMBase &fem);
            ~FEMPatchDirichletZero() override;

            PetscErrorCode ApplyBC(const Mat stiff_glb, Mat& stiff_inner, const Vec bc_glb, Vec& bc_inner) override;
            PetscErrorCode AugVec(const Vec vec_original, Vec& vec_aug) override;

        private:
            PetscErrorCode MakePatchNodeIndices() override;

    };
}


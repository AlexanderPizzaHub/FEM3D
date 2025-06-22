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

            const PetscInt GetPatchType() const {return patch_type_;};

        protected:
            femm::LagrangeFEMBase *fem_;

            PetscInt patch_type_, num_patch_nodes_;

            virtual PetscErrorCode MakePatchNodeIndices() = 0;
            
            PetscInt* patch_nodes_indices_; // 补丁节点的索引集

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
            PetscErrorCode AugVec(const Vec vec_original, Vec& vec_aug);

        private:

            IS inner_nodes_is_; // 内节点的索引集 
            PetscErrorCode MakePatchNodeIndices() override;

    };

    class FEMPatchDirichlet : public FEMPatchBase
    {
        public:
            FEMPatchDirichlet(femm::LagrangeFEMBase &fem, PetscInt patchlabel);
            ~FEMPatchDirichlet();

            PetscErrorCode BoundaryProject(PetscScalar (*func)(PetscScalar x, PetscScalar y, PetscScalar z), Vec &data);

            PetscErrorCode ApplyBC(const Mat stiff_glb, Mat& stiff_inner, const Vec bc_glb, Vec& bc_inner) override;
            PetscErrorCode AugVec(const Vec vec_original, Vec& vec_aug);


        private:
            PetscErrorCode MakePatchNodeIndices(PetscInt patchlabel);
            IS inner_nodes_is_;
            Vec data;
    };

    class FEMPatchNeumann : public FEMPatchBase 
    {
        public:
            FEMPatchNeumann(femm::LagrangeFEMBase &fem, PetscInt patchlabel);
            ~FEMPatchNeumann();

            PetscErrorCode ApplyBC(const Mat stiff_glb, Mat& stiff_inner, const Vec bc_glb, Vec& bc_inner) override;
            PetscErrorCode BoundaryProject(PetscScalar (*func)(PetscScalar x, PetscScalar y, PetscScalar z), Vec &data);
         
        private:
            PetscErrorCode MakePatchNodeIndices() override;
            Vec data;
    };
}


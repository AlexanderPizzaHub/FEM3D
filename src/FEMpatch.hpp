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
        FEMPatchBase(femm::LagrangeFEMBase &fem, PetscInt patchlabel);
        virtual ~FEMPatchBase() = default;

        PetscErrorCode BoundaryProject(PetscScalar (*func)(PetscScalar x, PetscScalar y, PetscScalar z));

        const PetscInt GetPatchType() const { return patch_type_; };
        const Vec GetData() const {return data;};

    protected:
        femm::LagrangeFEMBase *fem_;

        PetscInt patch_type_, num_patch_nodes_;
        PetscInt *patch_node_indices_; // 补丁节点的索引集
        Vec data;

        PetscErrorCode MakePatchNodeIndices(PetscInt patchlabel);
    };

    class FEMPatchDirichlet : public FEMPatchBase
    {
    public:
        FEMPatchDirichlet(femm::LagrangeFEMBase &fem, PetscInt patchlabel);
        ~FEMPatchDirichlet();

        PetscErrorCode ApplyBC(const Mat stiff_glb, Mat &stiff_inner, const Vec bc_glb, Vec &bc_inner);
        PetscErrorCode AugVec(const Vec vec_original, Vec &vec_aug);

    private:
        PetscInt num_inner_nodes_;
        PetscInt *inner_node_indices_;

        PetscErrorCode MakeInnerNodeIndices();
    };

    class FEMPatchNeumann : public FEMPatchBase
    {
    public:
        FEMPatchNeumann(femm::LagrangeFEMBase &fem, PetscInt patchlabel);
        ~FEMPatchNeumann();

        PetscErrorCode ApplyBC(Vec& rhs);
        PetscErrorCode ApplyBCBySource(Vec& source);

    };
}

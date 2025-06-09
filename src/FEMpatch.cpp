#include "FEMpatch.hpp"

namespace fempatch
{
    #pragma region base
    FEMPatchBase::FEMPatchBase(femm::LagrangeFEMBase &fem)
        : fem_(&fem) {
        // 构造函数
    }
    #pragma endregion



    #pragma region dirichletzero
    FEMPatchDirichletZero::FEMPatchDirichletZero(femm::LagrangeFEMBase &fem)
        : FEMPatchBase(fem) {
        // 构造函数
        MakePatchNodeIndices();
        std::cout << "FEMPatchDirichletZero initialized!" << std::endl;
    }

    FEMPatchDirichletZero::~FEMPatchDirichletZero() {
        // 析构函数
        //ISDestroy(&inner_nodes_is_);
        //ISDestroy(&patch_nodes_is_);
    }

    PetscErrorCode FEMPatchDirichletZero::ApplyBC(const Mat stiff_glb, Mat& stiff_inner, const Vec bc_glb, Vec& bc_inner) {
        // 应用Dirichlet边界条件
        PetscInt inner_size, patch_size;

        // 暂时处理makepatch中的bug        
        PetscCall(ISGetSize(inner_nodes_is_, &inner_size));
        PetscCall(ISGetSize(patch_nodes_is_, &patch_size));

        PetscCall(utils::MatSetup(inner_size, inner_size, stiff_inner));
        PetscCall(utils::VecSetup(inner_size, bc_inner));

        PetscCall(MatCreateSubMatrix(stiff_glb, inner_nodes_is_, inner_nodes_is_,MAT_INITIAL_MATRIX, &stiff_inner));
        PetscCall(VecGetSubVector(bc_glb, inner_nodes_is_, &bc_inner));
        
        return 0;
    }

    PetscErrorCode FEMPatchDirichletZero::AugVec(const Vec vec_original, Vec& vec_aug)
    {
        // 这个函数不该有的，临时添加
        const PetscInt *patchpts, *innerpts;
        PetscInt numpatchpts, numinnerpts;
        PetscCall(ISGetIndices(patch_nodes_is_, &patchpts));
        PetscCall(ISGetIndices(inner_nodes_is_, &innerpts));
        
        PetscCall(ISGetLocalSize(patch_nodes_is_, &numpatchpts));
        PetscCall(ISGetLocalSize(inner_nodes_is_, &numinnerpts));
        
        PetscScalar *inner_vals;
        PetscCall(VecGetArray(vec_original,&inner_vals));

        PetscCall(VecSetValues(vec_aug, numinnerpts,innerpts, inner_vals, INSERT_VALUES));
        return 0;
    }

    PetscErrorCode FEMPatchDirichletZero::MakePatchNodeIndices() {
        // 创建内节点和补丁节点的索引集
        const mesh::MeshDMPlex mesh = (fem_->GetMesh());
        inner_nodes_is_ = mesh.GetInnerIS();
        patch_nodes_is_ = mesh.GetBdryIS();
        //std::cout << inner_nodes_is_ << std::endl;
        return 0;
    }

    #pragma endregion

}
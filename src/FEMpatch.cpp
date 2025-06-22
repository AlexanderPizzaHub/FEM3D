#include "FEMpatch.hpp"

namespace fempatch
{
#pragma region base
    FEMPatchBase::FEMPatchBase(femm::LagrangeFEMBase &fem, PetscInt patchlabel)
        : fem_(&fem)
    {
        MakePatchNodeIndices(patchlabel);
        // 构造函数
    }

    PetscErrorCode FEMPatchBase::MakePatchNodeIndices(PetscInt patchlabel)
    {
        //DMLabel node_labels = fem_->GetNodeLabels();
        //PetscCall(DMLabelGetStratumIS(node_labels, patchlabel, &patch_nodes_is_));
        auto node_label_indices = fem_->GetNodeLabels();
        auto node_label_index = node_label_indices[patchlabel];
        num_patch_nodes_ = node_label_index.size();
        patch_node_indices_ = new PetscInt[num_patch_nodes_];
        for(PetscInt node_idx= 0 ; node_idx < num_patch_nodes_; node_idx++)
        {
            
            patch_node_indices_[node_idx] = node_label_index[node_idx];
        }
        return 0;
    }

    PetscErrorCode FEMPatchBase::BoundaryProject(PetscScalar (*func)(PetscScalar x, PetscScalar y, PetscScalar z))
    {
        PetscScalar pointdata;
        PetscScalar* coordinates;
        PetscInt node_idx;
        Vec node_coords = fem_->GetNodeCoords();
        std::cout << "$$$" <<std::endl;
        PetscCall(VecGetArray(node_coords,&coordinates));
        std::cout << "in bdproj:" << num_patch_nodes_<<std::endl;
        PetscCall(utils::VecSetup(num_patch_nodes_,data));
        //int s = patch_node_indices_.size();
        //std::cout << s<<": size"<<std::endl;
        for(PetscInt idx=0; idx<num_patch_nodes_; idx++)
        {
            //std::cout <<"@"<< idx <<std::endl;
            node_idx = patch_node_indices_[idx];
            //std::cout <<"@"<< node_idx <<std::endl;
            pointdata = func(
                coordinates[3*node_idx],
                coordinates[3*node_idx+1],
                coordinates[3*node_idx+2]
            );
            //std::cout <<"@"<< idx <<std::endl;
            PetscCall(VecSetValue(data, idx, pointdata, INSERT_VALUES));
            //std::cout <<"@"<< idx <<std::endl;
        }
        
        PetscCall(VecAssemblyBegin(data));
        PetscCall(VecAssemblyEnd(data));
        return 0;
    }
#pragma endregion

#pragma region dirichlet
    FEMPatchDirichlet::FEMPatchDirichlet(femm::LagrangeFEMBase &fem, PetscInt patchlabel)
        : FEMPatchBase(fem, patchlabel)
    {
        patch_type_ = 1;
        MakeInnerNodeIndices();
    }
    FEMPatchDirichlet::~FEMPatchDirichlet() {};

    PetscErrorCode FEMPatchDirichlet::MakeInnerNodeIndices()
    {
        std::cout << "in dir"<<std::endl;
        auto node_label_indices = fem_->GetNodeLabels();
        auto inner_node_label_index = node_label_indices[0]; // inner domain label
        std::cout << "111" <<std::endl;
    
        num_inner_nodes_ = inner_node_label_index.size();
        std::cout << "222:"<<num_inner_nodes_ <<std::endl;
        inner_node_indices_ = new PetscInt[num_inner_nodes_];
        for(PetscInt node_idx= 0 ; node_idx < num_inner_nodes_; node_idx++)
        {
            //std::cout <<"!!";
            inner_node_indices_[node_idx] = inner_node_label_index[node_idx];
        }
        //std::cout << "@@@" << std::endl;
        return 0;
    }


    PetscErrorCode FEMPatchDirichlet::ApplyBC(const Mat stiff_glb, Mat &stiff_inner, const Vec bc_glb, Vec &bc_inner)
    {
        // 应用Dirichlet边界条件
        /*
        
        还需要加入右端项修正
        
        */
        PetscCall(utils::MatSetup(num_inner_nodes_, num_inner_nodes_, stiff_inner));
        PetscCall(utils::VecSetup(num_inner_nodes_, bc_inner));

        IS inner_is;
        PetscCall(ISCreateGeneral(PETSC_COMM_WORLD,num_inner_nodes_,inner_node_indices_,PETSC_USE_POINTER,&inner_is));

        PetscCall(MatCreateSubMatrix(stiff_glb, inner_is, inner_is, MAT_INITIAL_MATRIX, &stiff_inner));
        PetscCall(VecGetSubVector(bc_glb, inner_is, &bc_inner));

        return 0;
    }

    PetscErrorCode FEMPatchDirichlet::AugVec(const Vec vec_original, Vec &vec_aug)
    {
        // 这个函数不该有的，临时添加
        PetscScalar *inner_vals, *patch_vals;
        PetscCall(VecGetArray(vec_original, &inner_vals));
        PetscCall(VecGetArray(data, &patch_vals));

        PetscCall(VecSetValues(vec_aug, num_inner_nodes_, inner_node_indices_, inner_vals, INSERT_VALUES));
        PetscCall(VecSetValues(vec_aug, num_patch_nodes_,patch_node_indices_, patch_vals,INSERT_VALUES));
        return 0;
    }

#pragma endregion

#pragma region neumann
FEMPatchNeumann::FEMPatchNeumann(femm::LagrangeFEMBase &fem, PetscInt patchlabel)
        : FEMPatchBase(fem, patchlabel)
    {
        patch_type_ = 2;
    }
FEMPatchNeumann::~FEMPatchNeumann() {};

PetscErrorCode FEMPatchNeumann::ApplyBC(Vec& rhs)
    {
        // 没写完，需要用到边界
        PetscScalar *patch_vals;
        PetscCall(VecGetArray(data, &patch_vals));
        
        PetscCall(VecSetValues(rhs, num_patch_nodes_, patch_node_indices_, patch_vals, ADD_VALUES));

        return 0;
    }

    PetscErrorCode FEMPatchNeumann::ApplyBCBySource(Vec& source)
    {
        // 没写完，需要用到边界
        PetscScalar *patch_vals;
        PetscCall(VecGetArray(data, &patch_vals));
        std::cout << num_patch_nodes_<<"!@#!"<<std::endl;
        PetscCall(VecSetValues(source, num_patch_nodes_, patch_node_indices_, patch_vals, ADD_VALUES));

        return 0;
    }

#pragma endregion

}
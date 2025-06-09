#include "FEMmachine.hpp"
#include <iostream>

namespace femm
{
#pragma region base
    LagrangeFEMBase::LagrangeFEMBase(mesh::MeshDMPlex &mesh, int order)
    {
        // 由于其他函数是纯虚函数，不能调用
        std::cout << "in base constructor!" << std::endl;
    }
    LagrangeFEMBase::~LagrangeFEMBase()
    {
        // 释放资源
        VecDestroy(&vert_indices_);
        VecDestroy(&vert_coords_);
        VecDestroy(&node_indices_);
        VecDestroy(&node_labels_);
        VecDestroy(&node_coords_);
        VecDestroy(&elem_indices_);
    }
#pragma endregion

#pragma region P1
    LagrangeP1FEM::LagrangeP1FEM(mesh::MeshDMPlex &mesh)
        : LagrangeFEMBase(mesh, 1) // P1有限元，阶数为1
    {
        mesh_ = &mesh;
        ExtractMeshData();
        MakeVert2NodeIndices();
        MakeVert2NodeCoords();
        MakeElemIndices();
        MakeElem2NodeMap();
        established = true; // 标记为已建立
    }
    LagrangeP1FEM::~LagrangeP1FEM()
    {
        // 释放资源
        if (established)
        {
            VecDestroy(&node_indices_);
            VecDestroy(&node_coords_);
            MatDestroy(&stiff_);
            MatDestroy(&rhs_);
        }
    }

    PetscErrorCode LagrangeP1FEM::ExtractMeshData()
    {
        PetscInt node_start, node_end;
        PetscScalar *global_coords, *coords;
        Vec coordinates;
        DM dm = mesh_->GetDM();

        PetscCall(DMPlexGetDepthStratum(dm, 0, &node_start, &node_end));
        PetscCall(DMGetCoordinatesLocal(dm, &coordinates));

        PetscCall(VecDuplicate(coordinates, &vert_coords_));
        PetscCall(VecCopy(coordinates, vert_coords_));

        for (PetscInt p = node_start; p < node_end; ++p)
        {
            VecSetValue(vert_indices_, p - node_start, p, INSERT_VALUES);
        }

        return 0; // 返回错误码，0表示成功
    }

    PetscErrorCode LagrangeP1FEM::MakeVert2NodeIndices()
    {
        // P1元节点与网格点相同,可以直接复制过去
        //PetscCall(VecDuplicate(vert_indices_, &node_indices_));
        //PetscCall(VecCopy(vert_indices_, node_indices_));
        node_indices_ = vert_indices_;
        PetscCall(VecGetSize(node_indices_, &num_nodes_));
        return 0; // 返回错误码，0表示成功
    }

    PetscErrorCode LagrangeP1FEM::MakeVert2NodeCoords()
    {
        // P1元节点与网格点相同,可以直接复制过去
        //PetscCall(VecDuplicate(vert_coords_, &node_coords_array));
        //PetscCall(VecCopy(vert_coords_, node_coords_array));
        node_coords_ = vert_coords_;
        return 0; // 返回错误码，0表示成功
    };

    PetscErrorCode LagrangeP1FEM::MakeElemIndices() 
    {
        PetscInt elem_start, elem_end;
        DM dm = mesh_->GetDM();

        PetscCall(DMPlexGetDepthStratum(dm, dim-1, &elem_start, &elem_end));
        PetscCall(utils::VecSetup(elem_end - elem_start, elem_indices_));

        for (PetscInt p = elem_start; p < elem_end; ++p)
        {
            VecSetValue(elem_indices_, p - elem_start, p, INSERT_VALUES);
        }

        PetscCall(VecGetSize(elem_indices_, &num_elements_));
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::GetElem2NodeIdx(PetscInt elem_idx, PetscInt* elem2node_idx)
    {
        mesh_->GetCell2VertIdx(elem_idx, elem2node_idx);
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::MakeElem2NodeMap() {
        elem2node_map_ = mesh_->GetCell2VertMap();
        return 0;
    };
    

    PetscErrorCode LagrangeP1FEM::AssembleStiff() {
        utils::MatSetup(num_nodes_, num_nodes_, stiff_);

        std::vector<PetscInt> elem2node_idx; // 四面体的四个顶点
        PetscScalar* J = new PetscScalar[dim * dim]; // Jacobian矩阵, row-first
        PetscScalar detJ = 0.0;
        PetscScalar detJinv = 0.0;
        PetscScalar *adjJ_nablaL0 = new PetscScalar[dim];
        PetscScalar *adjJ_nablaL1 = new PetscScalar[dim];
        PetscScalar *adjJ_nablaL2 = new PetscScalar[dim];
        PetscScalar *adjJ_nablaL3 = new PetscScalar[dim];
        PetscScalar *sub_stiff = new PetscScalar[(dim+1) * (dim+1)];
        PetscScalar *node_coords_array;

        PetscCall(VecGetArray(node_coords_, &node_coords_array));

        PetscScalar ref_vol = 1.0 / 6.0; // 四面体的参考体积
        
        PetscScalar tmpval = 0.0; // 储存临时变量

        for (PetscInt elem_idx = 0 ; elem_idx < num_elements_; ++elem_idx)
        {
            elem2node_idx = elem2node_map_[elem_idx];
            J[0] = node_coords_array[elem2node_idx[1] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[1] = node_coords_array[elem2node_idx[1] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[2] = node_coords_array[elem2node_idx[1] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[3] = node_coords_array[elem2node_idx[2] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[4] = node_coords_array[elem2node_idx[2] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[5] = node_coords_array[elem2node_idx[2] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[6] = node_coords_array[elem2node_idx[3] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[7] = node_coords_array[elem2node_idx[3] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[8] = node_coords_array[elem2node_idx[3] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            
            PetscCall(numerical::determinant33(J, detJ));
            detJinv = 1.0 / detJ;

            adjJ_nablaL1[0] = J[4] * J[8] - J[5] * J[7];
            adjJ_nablaL1[1] = J[5] * J[6] - J[3] * J[8];
            adjJ_nablaL1[2] = J[3] * J[7] - J[4] * J[6];

            adjJ_nablaL2[0] = J[2] * J[7] - J[1] * J[8];
            adjJ_nablaL2[1] = J[0] * J[8] - J[2] * J[6];
            adjJ_nablaL2[2] = J[1] * J[6] - J[0] * J[7];

            adjJ_nablaL3[0] = J[1] * J[5] - J[2] * J[4];
            adjJ_nablaL3[1] = J[2] * J[3] - J[0] * J[5];
            adjJ_nablaL3[2] = J[0] * J[4] - J[1] * J[3];

            adjJ_nablaL0[0] = -adjJ_nablaL1[0] - adjJ_nablaL2[0] - adjJ_nablaL3[0];
            adjJ_nablaL0[1] = -adjJ_nablaL1[1] - adjJ_nablaL2[1] - adjJ_nablaL3[1];
            adjJ_nablaL0[2] = -adjJ_nablaL1[2] - adjJ_nablaL2[2] - adjJ_nablaL3[2];

            // 计算刚度矩阵的子块
            PetscCall(numerical::vecinner(adjJ_nablaL0, adjJ_nablaL0, dim, tmpval));
            sub_stiff[0] = ref_vol * detJinv * tmpval;
            PetscCall(numerical::vecinner(adjJ_nablaL0, adjJ_nablaL1, dim, tmpval));
            sub_stiff[1] = ref_vol * detJinv * tmpval;
            sub_stiff[4] = sub_stiff[1]; 
            PetscCall(numerical::vecinner(adjJ_nablaL0, adjJ_nablaL2, dim, tmpval));
            sub_stiff[2] = ref_vol * detJinv * tmpval;
            sub_stiff[8] = sub_stiff[2];
            PetscCall(numerical::vecinner(adjJ_nablaL0, adjJ_nablaL3, dim, tmpval));
            sub_stiff[3] = ref_vol * detJinv * tmpval;
            sub_stiff[12] = sub_stiff[3];
            PetscCall(numerical::vecinner(adjJ_nablaL1, adjJ_nablaL1, dim, tmpval));
            sub_stiff[5] = ref_vol * detJinv * tmpval;
            PetscCall(numerical::vecinner(adjJ_nablaL1, adjJ_nablaL2, dim, tmpval));
            sub_stiff[6] = ref_vol * detJinv * tmpval;
            sub_stiff[9] = sub_stiff[6];
            PetscCall(numerical::vecinner(adjJ_nablaL1, adjJ_nablaL3, dim, tmpval));
            sub_stiff[7] = ref_vol * detJinv * tmpval;
            sub_stiff[13] = sub_stiff[7];
            PetscCall(numerical::vecinner(adjJ_nablaL2, adjJ_nablaL2, dim, tmpval));
            sub_stiff[10] = ref_vol * detJinv * tmpval;
            PetscCall(numerical::vecinner(adjJ_nablaL2, adjJ_nablaL3, dim, tmpval));
            sub_stiff[11] = ref_vol * detJinv * tmpval;
            sub_stiff[14] = sub_stiff[11];
            PetscCall(numerical::vecinner(adjJ_nablaL3, adjJ_nablaL3, dim, tmpval));
            sub_stiff[15] = ref_vol * detJinv * tmpval;

            const PetscInt *elem2node = elem2node_idx.data();
            PetscCall(MatSetValues(stiff_, dim+1, elem2node, dim+1, elem2node, sub_stiff, ADD_VALUES));
        }
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::AssembleRHS() {return 0;};

#pragma endregion
} // namespace femm
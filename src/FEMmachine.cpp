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
        MakeVert2NodeIndices();
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
        vert_indices_ = node_indices_;
        return 0; // 返回错误码，0表示成功
    }

    PetscErrorCode LagrangeP1FEM::MakeVert2NodeCoords()
    {
        // P1元节点与网格点相同,可以直接复制过去
        //PetscCall(VecDuplicate(vert_coords_, &node_coords_));
        //PetscCall(VecCopy(vert_coords_, node_coords_));
        vert_coords_ = node_coords_;
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
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::GetElem2NodeIdx(PetscInt elem_idx, PetscInt* elem2node_idx)
    {
        // 重点检查
        /*
        const PetscInt *edgecone, *nodecone;
        PetscInt node_start, node_end, edge_conesize, node_conesize;
        PetscInt neighbor_count;
        PetscInt *node_neighbors_tmp;

        PetscCall(DMPlexGetDepthStratum(dm, 0, &node_start, &node_end));
        PetscCall(DMPlexGetCone(dm, p, &edgecone));
        PetscCall(DMPlexGetConeSize(dm, p, &edge_conesize));

        assert(edge_conesize == 3);

        node_neighbors_tmp = new PetscInt[6];
        neighbor_count = 0;
        for (PetscInt i = 0; i < edge_conesize; i++)
        {
            PetscInt edge = edgecone[i];
            PetscCall(DMPlexGetConeSize(dm, edge, &node_conesize));
            PetscCall(DMPlexGetCone(dm, edge, &nodecone));

            assert(node_conesize == 2);
            for (PetscInt j = 0; j < node_conesize; j++)
            {
                PetscInt node = nodecone[j];
                assert(node >= node_start && node < node_end);
                node_neighbors_tmp[neighbor_count] = node;
                neighbor_count++;
            }
        }
        assert(neighbor_count == 6);
        PetscCall(PetscSortInt(neighbor_count, node_neighbors_tmp));

        for (PetscInt i = 0; i < 3; i++)
        {
            node_neighbors[i] = node_neighbors_tmp[2 * i];
        }
            */
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::MakeElem2VertMap()
    {
        PetscInt elem_start, elem_end;
        PetscInt *elem2vert_idx;

        DM dm = mesh_->GetDM();
        PetscCall(DMPlexGetDepthStratum(dm, dim-1, &elem_start, &elem_end));
        elem2vert_map_.reserve(elem_end-elem_start);
        elem2vert_idx = new PetscInt[4]; // 四面体的四个顶点

        for (PetscInt p_elem = elem_start; p_elem < elem_end; ++p_elem)
        {
            PetscCall(GetElem2NodeIdx(p_elem,elem2vert_idx));
            std::vector<PetscInt> vert_indices(4);
            for (PetscInt i = 0; i < 4; i++)
            {
                vert_indices[i] = elem2vert_idx[i];
            }
            elem2vert_map_.push_back(vert_indices);
        }
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::MakeElem2NodeMap() {
        

        // 我们建elem2


    };
    

    PetscErrorCode LagrangeP1FEM::AssembleStiff() {};

    PetscErrorCode LagrangeP1FEM::AssembleRHS() {};

#pragma endregion
} // namespace femm
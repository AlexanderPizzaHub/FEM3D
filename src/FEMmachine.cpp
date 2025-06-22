#include "FEMmachine.hpp"
#include <iostream>

namespace femm
{
#pragma region base
    LagrangeFEMBase::LagrangeFEMBase(mesh::MeshDMPlex &mesh, int order)
    {
        // 由于其他函数是纯虚函数，不能调用
        order_ = order;
        mesh_ = &mesh;
        std::cout << "in base constructor!" << std::endl;
    }
    LagrangeFEMBase::~LagrangeFEMBase()
    {

        if (established)
        {
            VecDestroy(&vert_indices_);
            VecDestroy(&vert_coords_);
            VecDestroy(&node_indices_);
            VecDestroy(&node_coords_);
            VecDestroy(&elem_indices_);
        }
        // 释放资源
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
        std::cout << "FEM initialization done!" << std::endl;
        // established = true; // 标记为已建立
    }
    LagrangeP1FEM::~LagrangeP1FEM()
    {
        // 释放资源
        if (established)
        {
            VecDestroy(&node_indices_);
            VecDestroy(&node_coords_);
            MatDestroy(&stiff_);
            VecDestroy(&rhs_);
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

        PetscCall(utils::VecSetup(node_end - node_start, vert_indices_));

        for (PetscInt p = node_start; p < node_end; ++p)
        {
            VecSetValue(vert_indices_, p - node_start, p, INSERT_VALUES);
        }

        DMLabel vertex_label = mesh_->GetLabel();
        PetscInt num_patches;
        PetscCall(DMLabelGetNumValues(vertex_label, &num_patches));

        IS patchIS;
        PetscInt patch_size;
        const PetscInt *vertex_indices;
        for(PetscInt patch_idx = 0; patch_idx < num_patches; patch_idx++)
        {
            PetscCall(DMLabelGetStratumIS(vertex_label, patch_idx, &patchIS));

            PetscCall(ISGetLocalSize(patchIS, &patch_size));
            //PetscCall(ISShift(patchIS, -node_start, patchIS));

            PetscCall(ISGetIndices(patchIS, &vertex_indices));
    
            std::vector<PetscInt> node_label_idx;
            for(PetscInt subidx=0; subidx<patch_size;subidx++)
            {
                node_label_idx.push_back(vertex_indices[subidx]);
            }
            
            node_label_indices_.push_back(node_label_idx);
        }
        //std::cout << num_patches<<std::endl;
        return 0; // 返回错误码，0表示成功
    }

    PetscErrorCode LagrangeP1FEM::MakeVert2NodeIndices()
    {
        // P1元节点与网格点相同,可以直接复制过去
        // PetscCall(VecDuplicate(vert_indices_, &node_indices_));
        // PetscCall(VecCopy(vert_indices_, node_indices_));
        node_indices_ = vert_indices_;
        PetscCall(VecGetSize(node_indices_, &num_nodes_));
        return 0; // 返回错误码，0表示成功
    }

    PetscErrorCode LagrangeP1FEM::MakeVert2NodeCoords()
    {
        // P1元节点与网格点相同,可以直接复制过去
        // PetscCall(VecDuplicate(vert_coords_, &node_coords_array));
        // PetscCall(VecCopy(vert_coords_, node_coords_array));
        node_coords_ = vert_coords_;
        return 0; // 返回错误码，0表示成功
    };

    PetscErrorCode LagrangeP1FEM::MakeElemIndices()
    {
        PetscInt elem_start, elem_end;
        DM dm = mesh_->GetDM();

        PetscCall(DMPlexGetDepthStratum(dm, dim, &elem_start, &elem_end));
        PetscCall(utils::VecSetup(elem_end - elem_start, elem_indices_));

        for (PetscInt p = elem_start; p < elem_end; ++p)
        {
            VecSetValue(elem_indices_, p - elem_start, p, INSERT_VALUES);
        }

        PetscCall(VecGetSize(elem_indices_, &num_elements_));
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::GetElem2NodeIdx(PetscInt elem_idx, PetscInt *elem2node_idx)
    {
        mesh_->GetCell2VertIdx(elem_idx, elem2node_idx);
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::MakeElem2NodeMap()
    {
        PetscInt node_start, node_end;
        DM dm = mesh_->GetDM();
        PetscCall(DMPlexGetDepthStratum(dm, 0, &node_start, &node_end));

        elem2node_map_ = mesh_->GetCell2VertMap();
        for (PetscInt i = 0; i < elem2node_map_.size(); ++i)
        {
            for (PetscInt j = 0; j < elem2node_map_[i].size(); ++j)
            {
                // 将网格点索引转换为节点索引
                elem2node_map_[i][j] -= node_start;
            }
        }
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::AssembleStiff()
    {
        utils::MatSetup(num_nodes_, num_nodes_, stiff_);

        std::vector<PetscInt> elem2node_idx;         // 四面体的四个顶点
        PetscScalar *J = new PetscScalar[dim * dim]; // Jacobian矩阵, row-first
        PetscScalar detJ = 0.0;
        PetscScalar detJinv = 0.0;
        PetscScalar *adjJ_nablaL0 = new PetscScalar[dim];
        PetscScalar *adjJ_nablaL1 = new PetscScalar[dim];
        PetscScalar *adjJ_nablaL2 = new PetscScalar[dim];
        PetscScalar *adjJ_nablaL3 = new PetscScalar[dim];
        PetscScalar *sub_stiff = new PetscScalar[(dim + 1) * (dim + 1)];
        PetscScalar *node_coords_array;

        PetscCall(VecGetArray(node_coords_, &node_coords_array));

        PetscScalar ref_vol = 1.0 / 6.0; // 四面体的参考体积

        PetscScalar tmpval = 0.0; // 储存临时变量
        std::cout << "assembling stiffness matrix..." << std::endl;
        std::cout << "num_elements_: " << num_elements_ << std::endl;

        for (PetscInt elem_idx = 0; elem_idx < num_elements_; ++elem_idx)
        {
            elem2node_idx = elem2node_map_[elem_idx];
            
            /*
            J[0] = node_coords_array[elem2node_idx[1] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[1] = node_coords_array[elem2node_idx[1] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[2] = node_coords_array[elem2node_idx[1] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[3] = node_coords_array[elem2node_idx[2] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[4] = node_coords_array[elem2node_idx[2] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[5] = node_coords_array[elem2node_idx[2] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[6] = node_coords_array[elem2node_idx[3] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[7] = node_coords_array[elem2node_idx[3] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[8] = node_coords_array[elem2node_idx[3] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            */
            

           J[0] = node_coords_array[elem2node_idx[1] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[1] = node_coords_array[elem2node_idx[2] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[2] = node_coords_array[elem2node_idx[3] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[3] = node_coords_array[elem2node_idx[1] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[4] = node_coords_array[elem2node_idx[2] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[5] = node_coords_array[elem2node_idx[3] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[6] = node_coords_array[elem2node_idx[1] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[7] = node_coords_array[elem2node_idx[2] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[8] = node_coords_array[elem2node_idx[3] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            
            
            numerical::determinant33(J, detJ);
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
            numerical::vecinner(adjJ_nablaL0, adjJ_nablaL0, dim, tmpval);
            sub_stiff[0] = ref_vol * detJinv * tmpval;
            if (sub_stiff[0] < 0.0)
            {
                // 说明标号顺序导致对角出负。处理网格麻烦，所以我们在这里手动处理
                sub_stiff[0] *= -1.0;
                detJinv *= -1.0;
            }
            numerical::vecinner(adjJ_nablaL0, adjJ_nablaL1, dim, tmpval);
            sub_stiff[1] = ref_vol * detJinv * tmpval;
            sub_stiff[4] = sub_stiff[1];
            numerical::vecinner(adjJ_nablaL0, adjJ_nablaL2, dim, tmpval);
            sub_stiff[2] = ref_vol * detJinv * tmpval;
            sub_stiff[8] = sub_stiff[2];
            numerical::vecinner(adjJ_nablaL0, adjJ_nablaL3, dim, tmpval);
            sub_stiff[3] = ref_vol * detJinv * tmpval;
            sub_stiff[12] = sub_stiff[3];
            numerical::vecinner(adjJ_nablaL1, adjJ_nablaL1, dim, tmpval);
            sub_stiff[5] = ref_vol * detJinv * tmpval;
            numerical::vecinner(adjJ_nablaL1, adjJ_nablaL2, dim, tmpval);
            sub_stiff[6] = ref_vol * detJinv * tmpval;
            sub_stiff[9] = sub_stiff[6];
            numerical::vecinner(adjJ_nablaL1, adjJ_nablaL3, dim, tmpval);
            sub_stiff[7] = ref_vol * detJinv * tmpval;
            sub_stiff[13] = sub_stiff[7];
            numerical::vecinner(adjJ_nablaL2, adjJ_nablaL2, dim, tmpval);
            sub_stiff[10] = ref_vol * detJinv * tmpval;
            numerical::vecinner(adjJ_nablaL2, adjJ_nablaL3, dim, tmpval);
            sub_stiff[11] = ref_vol * detJinv * tmpval;
            sub_stiff[14] = sub_stiff[11];
            numerical::vecinner(adjJ_nablaL3, adjJ_nablaL3, dim, tmpval);
            sub_stiff[15] = ref_vol * detJinv * tmpval;

            const PetscInt *elem2node = elem2node_idx.data();
            // std::cout << "elem_idx: " << elem_idx << ", elem2node: ";
            /*
            std::cout << "local stiff matrix: " << std::endl;
            for (PetscInt i = 0; i < dim + 1; ++i)
            {
                for (PetscInt j = 0; j < dim + 1; ++j)
                {
                    std::cout << sub_stiff[i * (dim + 1) + j] << " ";
                }
                std::cout << std::endl;
            }
                */
            /*
            std::cout << "element: " << elem_idx<<std::endl;
            std::cout << "substiff:" << std::endl;
            for(int i =0;i<4;i++)
            {
                for(int j =0;j<4;j++)
                {
                    std::cout << sub_stiff[4*i+j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "coords: "<< std::endl;
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<3;j++)
                {
                    std::cout << node_coords_array[elem2node_idx[i] * dim+j]<<" ";
                }
                std::cout << std::endl;
            }
            std::cout << "J adjoint: "<<std::endl;
            std::cout << adjJ_nablaL1[0] << " " << adjJ_nablaL1[1] << " "<<adjJ_nablaL1[2] << std::endl;
            std::cout << adjJ_nablaL2[0] << " " << adjJ_nablaL2[1] << " "<<adjJ_nablaL2[2] << std::endl;
            std::cout << adjJ_nablaL3[0] << " " << adjJ_nablaL3[1] << " "<<adjJ_nablaL3[2] << std::endl;
            */

            PetscCall(MatSetValues(stiff_, dim + 1, elem2node, dim + 1, elem2node, sub_stiff, ADD_VALUES));
        }
        // std::cout << "AssembleStiff done!" << std::endl;
        PetscCall(MatAssemblyBegin(stiff_, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff_, MAT_FINAL_ASSEMBLY));


        return 0;
    };

    PetscErrorCode LagrangeP1FEM::AssembleStiffDeg()
    {
        utils::MatSetup(num_nodes_, num_nodes_, stiff_);

        std::vector<PetscInt> elem2node_idx; // 四面体的四个顶点

        PetscScalar *sub_stiff = new PetscScalar[(dim + 1) * (dim + 1)];

        std::cout << "assembling stiffness deg matrix..." << std::endl;
        std::cout << "num_elements_: " << num_elements_ << std::endl;

        for (PetscInt elem_idx = 0; elem_idx < num_elements_; ++elem_idx)
        {
            elem2node_idx = elem2node_map_[elem_idx];

            sub_stiff[0] = 3;
            sub_stiff[1] = -1;
            sub_stiff[4] = -1;

            sub_stiff[2] = -1;
            sub_stiff[8] = -1;

            sub_stiff[3] = -1;
            sub_stiff[12] = -1;

            sub_stiff[5] = 3;

            sub_stiff[6] = -1;
            sub_stiff[9] = -1;

            sub_stiff[7] = -1;
            sub_stiff[13] = -1;

            sub_stiff[10] = 3;

            sub_stiff[11] = -1;
            sub_stiff[14] = -1;

            sub_stiff[15] = 3;

            const PetscInt *elem2node = elem2node_idx.data();
            // std::cout << "elem_idx: " << elem_idx << ", elem2node: ";

            PetscCall(MatSetValues(stiff_, dim + 1, elem2node, dim + 1, elem2node, sub_stiff, ADD_VALUES));
        }
        // std::cout << "AssembleStiff done!" << std::endl;
        PetscCall(MatAssemblyBegin(stiff_, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff_, MAT_FINAL_ASSEMBLY));
        return 0;
    };



    PetscErrorCode LagrangeP1FEM::AssembleRHS(Vec data)
    {
        PetscScalar source, detJ;
        PetscScalar *data_section, *add_values, *J, *node_coords_array;
        PetscInt *elem2node_idx;

        data_section = new PetscScalar[4];
        add_values = new PetscScalar[4];
        J = new PetscScalar[9];

        PetscCall(VecGetArray(node_coords_, &node_coords_array));
        PetscCall(utils::VecSetup(num_nodes_,rhs_));
        PetscCall(VecSet(rhs_, 0.0));
        for (PetscInt p = 0; p < num_elements_; p++)
        {
            elem2node_idx = elem2node_map_[p].data();

            J[0] = node_coords_array[elem2node_idx[1] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[1] = node_coords_array[elem2node_idx[1] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[2] = node_coords_array[elem2node_idx[1] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[3] = node_coords_array[elem2node_idx[2] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[4] = node_coords_array[elem2node_idx[2] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[5] = node_coords_array[elem2node_idx[2] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];
            J[6] = node_coords_array[elem2node_idx[3] * dim] - node_coords_array[elem2node_idx[0] * dim];
            J[7] = node_coords_array[elem2node_idx[3] * dim + 1] - node_coords_array[elem2node_idx[0] * dim + 1];
            J[8] = node_coords_array[elem2node_idx[3] * dim + 2] - node_coords_array[elem2node_idx[0] * dim + 2];

            numerical::determinant33(J, detJ);

            PetscCall(VecGetValues(data, 4, elem2node_idx, data_section));
            //source = (data_section[0] + data_section[1] + data_section[2] + data_section[3]) / 4.0;
            source = data_section[0];

            

            if(detJ < 0){detJ *= -1.0;};

            source *= (detJ / 24.0);
            //source *= (1.0/24.0/detJ);
            //std::cout << source << std::endl;

            add_values[0] = source;
            add_values[1] = source;
            add_values[2] = source;
            add_values[3] = source;
            PetscCall(VecSetValues(rhs_, 4, elem2node_idx, add_values, ADD_VALUES));
        }
        PetscCall(VecAssemblyBegin(rhs_));
        PetscCall(VecAssemblyEnd(rhs_));
        //PetscCall(VecView(rhs_,PETSC_VIEWER_STDOUT_WORLD));
        
        return 0;
    };

    PetscErrorCode LagrangeP1FEM::DomainProject(PetscScalar (*func)(PetscScalar x, PetscScalar y, PetscScalar z), Vec &data)
    {
        PetscScalar *node_coords_array;
        PetscScalar data_point;

        PetscCall(utils::VecSetup(num_nodes_, data));
        PetscCall(VecGetArray(node_coords_, &node_coords_array));

        for (PetscInt i = 0; i < num_nodes_; i++)
        {
            data_point = func(node_coords_array[dim * i],
                              node_coords_array[dim * i + 1],
                              node_coords_array[dim * i + 2]);
            PetscCall(VecSetValue(data,i,data_point,INSERT_VALUES));
        }
        PetscCall(VecAssemblyBegin(data));
        PetscCall(VecAssemblyEnd(data));
        return 0;
    };

#pragma endregion
} // namespace femm
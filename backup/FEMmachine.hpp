/*
此模块负责各类FEM的实现。在每一类FEM中，包含：
    1. 将其他接口读入的网格信息转化为FEM所需的网格数据结构。 
    2. 相关组矩阵操作
*/

/*
remark: 这是将femmesh和fem区分开来的版本
*/

#pragma once

#include <petsc.h>
#include <string>
#include <vector>
#include <array>

#include "mesh.hpp"

namespace femm
{
    #pragma region base
    class LagrangeMeshBase
    {
        /*
        所有的索引默认从0开始
        假设矩阵向量数据结构底座由petsc提供
        */

        public: 
            // 从DMPlex中构造FEM需要的数据结构。未来可能需要从其他格式中构造，重载构造函数即可
            LagrangeMeshBase(mesh::MeshDMPlex &mesh, int order);

            ~LagrangeMeshBase();

            // 将网格点(vertex)的索引转换为有限元节点(node)的索引。传参用vec还是指针需要斟酌
            //virtual PetscErrorCode MakeVert2NodeIndices( PetscInt* &vert_indices)  = 0;
            virtual PetscErrorCode MakeVert2NodeIndices()  = 0;

            // 将网格点(vertex)的坐标转换为有限元节点(node)的坐标
            //virtual PetscErrorCode MakeVert2NodeCoords( std::vector<PetscScalar> &vert_coords)  = 0;
            virtual PetscErrorCode MakeVert2NodeCoords()  = 0;

            // 将网格单元(cell)的索引转换为有限元单元的索引。通常是相同的
            //virtual PetscErrorCode MakeElemIndices( std::vector<PetscInt> &element_indices)  = 0;
            virtual PetscErrorCode MakeElemIndices()  = 0;

            virtual PetscErrorCode MakeElem2NodeMap()  = 0;
            virtual PetscErrorCode MakeElem2VertMap()  = 0;

        protected:
            PetscInt order_; // 有限元的阶数

            // 我们仍然需要用到网格点的坐标，用于积分计算
            Vec vert_indices_; // 网格点的索引
            
            Vec vert_coords_; // 网格点的坐标

            Vec node_indices_; // 节点索引
            Vec node_labels_; // 节点标签，用于标记节点所属的元素类型

            Vec node_coords_; // 节点坐标

            Vec elem_indices_; // 单元索引

            // 这个还不清楚用petsc提供的数据格式，如何实现更好
            std::vector<std::vector<PetscInt>> elem2node_map_; // 元素到节点的映射
            std::vector<std::vector<PetscInt>> elem2vert_map_; // 元素到网格点的映射

            PetscInt num_nodes_; // 节点的数量
            PetscInt num_elements_; // 元素的数量

        private:
            mesh::MeshDMPlex *mesh_; // 网格数据结构

    };

    class LagrangeFEMBase
    {
        /*
        矩阵均为组全域
        */
        public:
            LagrangeFEMBase(LagrangeMeshBase *femmesh, int order);

            virtual ~LagrangeFEMBase();

            // 组装全局矩阵
            virtual PetscErrorCode AssembleStiff() = 0;

            virtual PetscErrorCode AssembleRHS() = 0;

            // 获取有限元网格数据结构
            LagrangeMeshBase* GetMesh()  { return femmesh_; }

        private:
            LagrangeMeshBase *femmesh_; // FEM网格数据结构
            PetscInt order_; // 有限元的阶数
            
            Mat stiff_; // 全局刚度矩阵
            Mat rhs_; // 全局右端项矩阵
    };
    #pragma endregion

    #pragma region femmesh
    class LagrangeP1FEM : public 
    {

    }

    #pragma endregion
}
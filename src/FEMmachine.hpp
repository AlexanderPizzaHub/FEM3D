/*
此模块负责各类FEM的实现。在每一类FEM中，包含：
    1. 将其他接口读入的网格信息转化为FEM所需的网格数据结构。 
    2. 相关组矩阵操作
*/

#pragma once

#include <petsc.h>
#include <string>
#include <vector>
#include <array>

#include "mesh.hpp"
#include "numericaltools.hpp"
#include "utils.hpp"

namespace femm
{
    #pragma region base
    class LagrangeFEMBase
    {
        /*
        所有的索引默认从0开始
        假设矩阵向量数据结构底座由petsc提供
        矩阵均为组全域
        */
        public:
            LagrangeFEMBase(mesh::MeshDMPlex &mesh, int order);

            virtual ~LagrangeFEMBase();

            /*-------------
            Mesh部分
            ---------------*/
            // 将网格点序号，坐标等提取出来
            virtual PetscErrorCode ExtractMeshData() = 0;

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

            /*--------------
            FEM部分
            ----------------*/

            // 组装全局矩阵
            virtual PetscErrorCode AssembleRHS(Vec data) = 0;

            virtual PetscErrorCode DomainProject(PetscScalar (*func)(PetscScalar x, PetscScalar y, PetscScalar z), Vec &data) = 0;

            mesh::MeshDMPlex &GetMesh() const { return *mesh_; }
            PetscInt GetNumNodes() const {return num_nodes_;}
            const std::vector<std::vector<PetscInt>>  GetNodeLabels() const {return node_label_indices_;};
            

        protected:
            /*----
            Mesh部分
            ----*/
            mesh::MeshDMPlex *mesh_; // 网格数据结构

            PetscInt order_; // 有限元的阶数

            // 我们仍然需要用到网格点的坐标，用于积分计算
            Vec vert_indices_; // 网格点的索引
            
            Vec vert_coords_; // 网格点的坐标

            Vec node_indices_; // 节点索引
            std::vector<std::vector<PetscInt>> node_label_indices_;
            //DMLabel node_labels_;
            //Vec node_labels_; // 节点标签，用于标记节点所属的元素类型
            

            Vec node_coords_; // 节点坐标

            Vec elem_indices_; // 单元索引

            // 这个还不清楚用petsc提供的数据格式，如何实现更好
            std::vector<std::vector<PetscInt>> elem2node_map_; // 元素到节点的映射

            PetscInt num_nodes_; // 节点的数量
            PetscInt num_elements_; // 元素的数量

            bool established = false;

    };
    #pragma endregion


    #pragma region femmesh
    class LagrangeP1FEM : public LagrangeFEMBase
    {
        /*
        P1有限元，线性有限元
        */
        public:
            LagrangeP1FEM(mesh::MeshDMPlex &mesh);

            ~LagrangeP1FEM();

            // 提取网格数据
            PetscErrorCode ExtractMeshData() override;

            // 将网格点(vertex)的索引转换为有限元节点(node)的索引
            PetscErrorCode MakeVert2NodeIndices()  override;

            // 将网格点(vertex)的坐标转换为有限元节点(node)的坐标
            PetscErrorCode MakeVert2NodeCoords()  override;

            // 将网格单元(cell)的索引转换为有限元单元的索引
            PetscErrorCode MakeElemIndices()  override;

            PetscErrorCode MakeElem2NodeMap()  override;

            // 组装全局刚度矩阵
            PetscErrorCode AssembleStiff();
            PetscErrorCode AssembleStiffDeg();

            // 组装全局质量矩阵

            // 组装全局右端项矩阵
            PetscErrorCode AssembleRHS(Vec data) override;

            PetscErrorCode GetElem2NodeIdx(PetscInt elem_idx, PetscInt* elem2node_idx);

            PetscErrorCode DomainProject(PetscScalar (*func)(PetscScalar x, PetscScalar y,PetscScalar z), Vec &data) override;

            Mat GetStiff() const { return stiff_; }
            Vec GetRHS() const { return rhs_; }
            
            
        private:
            /*----
                FEM部分
                ----*/
                
                Mat stiff_; // 全局刚度矩阵
                Vec rhs_; // 全局右端项矩阵

    };
    #pragma endregion
}
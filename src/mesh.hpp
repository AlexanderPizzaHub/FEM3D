/*
此模块负责各种网格格式，例如.msh, .cgns文件的读入。读入的数据存储为一个给定的，由接口提供的中间数据结构。

我们采用PETSc提供的DMPlex作为中间数据结构。
*/

#pragma once

#include <petsc.h>
#include <string>
#include <vector>
#include <array>
#include "numericaltools.hpp"
#include "utils.hpp"

namespace mesh
{
    class MeshDMPlex
    {
    public:
        MeshDMPlex(const char filename[]);

        // 第三种读取方式，我们手动给DMPlex标记边界。网格边界应该由mesh文件指定，因此该功能计划在未来废弃。
        // MeshDMPlex(const char filename[], const std::string domain_type); // 暂时默认为cube

        ~MeshDMPlex();

        // 根据label获取对应的网格点序号。
        PetscErrorCode GetISbyLabel(PetscInt label_index, IS &is);

        const DM &GetDM() const { return dm_; }
        const DMLabel &GetLabel() const { return label_; }
        const IS &GetBdryIS() const { return bdryIS_; }
        const IS &GetInnerIS() const { return innerIS_; }

    private:
        // 用于手动标记边界的函数。计划在未来废弃。
        PetscErrorCode MarkCubeBoundary( PetscScalar left, PetscScalar right, PetscScalar front, PetscScalar back, PetscScalar top, PetscScalar bottom);

        DM dm_;
        DMLabel label_; // 网格点的标签，用于区分边界，多介质区等。
        IS bdryIS_, innerIS_; // 边界和内点的索引集，用于标记边界和内点。这里没有考虑到混合边界的情况。该功能计划在未来改进。

    };
}
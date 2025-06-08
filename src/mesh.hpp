/*
此模块负责各种网格格式，例如.msh, .cgns文件的读入。读入的数据存储为一个给定的，由接口提供的中间数据结构。

我们采用PETSc提供的DMPlex作为中间数据结构。
*/

#pragma once

#include <petsc.h>
#include <string>
#include <vector>
#include <array>

namespace mesh
{
    class MeshDMPlex
    {
    public:
        MeshDMPlex(DM dm);
        MeshDMPlex(const char filename[]);

        // 第三种读取方式，我们手动给DMPlex标记边界。网格边界应该由mesh文件指定，因此该功能计划在未来废弃。
        MeshDMPlex(const char filename[], const std::string domain_type);

        ~MeshDMPlex();

        // 用于手动标记边界的函数。计划在未来废弃。
        PetscErrorCode MarkRectBoundary(DM dm, PetscScalar left, PetscScalar right, PetscScalar top, PetscScalar bottom, DMLabel &label);

        PetscErrorCode MarkAnnulusBoundary(DM dm, PetscScalar inner_radius, PetscScalar outer_radius, DMLabel &label);

        // 根据label获取对应的网格点序号。
        PetscErrorCode GetISbyLabel(DMLabel label, IS &is);

        const DM &GetDM() const { return dm; }
        const DMLabel &GetLabel() const { return label; }
        const IS &GetBdryIS() const { return bdryIS; }
        const IS &GetInnerIS() const { return innerIS; }

    private:
        DM dm;
        DMLabel label; // 网格点的标签，用于区分边界，多介质区等。
        IS bdryIS, innerIS; // 边界和内点的索引集，用于标记边界和内点。这里没有考虑到混合边界的情况。该功能计划在未来改进。

    };
}
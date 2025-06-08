#include "test.hpp"

namespace test
{
    PetscErrorCode TestMeshDMPlex()
    {
        // 测试MeshDMPlex的功能
        const char *filename = "../meshfile/cube.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);

        // 测试获取DM和Label
        const DM &dm = mesh.GetDM();
        const DMLabel &label = mesh.GetLabel();

        // 测试获取边界和内点的索引集
        const IS &bdryIS = mesh.GetBdryIS();
        const IS &innerIS = mesh.GetInnerIS();

        // 输出一些信息以验证功能
        PetscInt bdrySize, innerSize;
        ISGetSize(bdryIS, &bdrySize);
        //innerSize = 0; // 初始化innerSize为0
        ISGetSize(innerIS, &innerSize);
        
        PetscPrintf(PETSC_COMM_WORLD, "Boundary size: %d, Inner size: %d\n", bdrySize, innerSize);

        return 0;
    }

    PetscErrorCode TestFEMMachine()
    {
        // 测试FEM机器的功能
        const char *filename = "../meshfile/cube.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);
        
        // 创建P1有限元实例
        femm::LagrangeP1FEM fem(mesh);
        
        PetscPrintf(PETSC_COMM_WORLD, "FEM machine test completed successfully.\n");

        return 0;
    }
}
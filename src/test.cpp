#include "test.hpp"

namespace test
{
    PetscErrorCode TestMeshDMPlex()
    {
        // 测试MeshDMPlex的功能
        const char *filename = "../meshfile/cubefine.msh"; // 测试文件
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

    PetscErrorCode TestMeshVertMap()
    {
        const char *filename = "../meshfile/cubefine.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);

        auto cell2vert_map = mesh.GetCell2VertMap();
        std::cout << "Cell to Vertex Map:" << std::endl;
        for (int i = 0; i < cell2vert_map.size(); ++i)
        {
            std::cout << "Cell " << i << ": ";
            for (const auto &vert : cell2vert_map[i])
            {
                std::cout << vert << " ";
            }
            std::cout << std::endl;
        }

        return 0;
    }

    PetscErrorCode TestFEMMachine()
    {
        // 测试FEM机器的功能
        const char *filename = "../meshfile/cubefine.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);
        
        // 创建P1有限元实例
        femm::LagrangeP1FEM fem(mesh);
        
        PetscPrintf(PETSC_COMM_WORLD, "FEM machine test completed successfully.\n");

        return 0;
    }

    PetscErrorCode TestStiff()
    {
        const char *filename = "../meshfile/cubefine.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);
        
        // 创建P1有限元实例
        femm::LagrangeP1FEM fem(mesh);
        
        fem.AssembleStiff();
        Mat stiff = fem.GetStiff();
        PetscCall(MatAssemblyBegin(stiff, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff, MAT_FINAL_ASSEMBLY));
        PetscCall(MatView(stiff, PETSC_VIEWER_STDOUT_WORLD));
        /*
        for (int i =0; i< 63; i++)
        {
            PetscScalar val;
            PetscCall(MatGetValue(stiff, i, i, &val));
            PetscPrintf(PETSC_COMM_WORLD, "Stiffness matrix value at (%d, %d): %f\n", i, i, val);
        }
            */
        return 0;
    }

    PetscErrorCode TestStiffDeg()
    {
        const char *filename = "../meshfile/cubefine.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);
        
        // 创建P1有限元实例
        femm::LagrangeP1FEM fem(mesh);
        
        fem.AssembleStiffDeg();
        Mat stiff = fem.GetStiff();
        PetscCall(MatAssemblyBegin(stiff, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff, MAT_FINAL_ASSEMBLY));
        //PetscCall(MatView(stiff, PETSC_VIEWER_STDOUT_WORLD));
        /*
        for (int i =0; i< 63; i++)
        {
            PetscScalar val;
            PetscCall(MatGetValue(stiff, i, i, &val));
            PetscPrintf(PETSC_COMM_WORLD, "Stiffness matrix value at (%d, %d): %f\n", i, i, val);
        }
            */
            
        return 0;
    }

    PetscErrorCode TestPatch()
    {
        const char *filename = "../meshfile/cubefine.msh"; // 测试文件
        mesh::MeshDMPlex mesh(filename);
        
        // 创建P1有限元实例
        femm::LagrangeP1FEM fem(mesh);
        
        fem.AssembleStiff();
        Mat stiff = fem.GetStiff();
        PetscCall(MatAssemblyBegin(stiff, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff, MAT_FINAL_ASSEMBLY));

        fempatch::FEMPatchDirichletZero patch(fem);

        Vec bc_glb;
        PetscCall(utils::VecSetup(fem.GetNumNodes(), bc_glb));

        Mat stiff_inner;
        Vec bc_inner;

        patch.ApplyBC(stiff, stiff_inner, bc_glb, bc_inner);
        PetscCall(MatAssemblyBegin(stiff_inner, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff_inner, MAT_FINAL_ASSEMBLY));
        PetscCall(MatView(stiff_inner, PETSC_VIEWER_STDOUT_WORLD));
    
        return 0;
    }

    PetscErrorCode TestSolver()
    {
        const char *filename = "../meshfile/cubefine1.msh"; // 测试文件
        std::__1::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

        mesh::MeshDMPlex mesh(filename);

        std::__1::chrono::steady_clock::time_point now1 = std::chrono::steady_clock::now();
        double durat1 = std::chrono::duration_cast<std::chrono::milliseconds>(now1-start).count();
        std::cout << "mesh loaded! time: "<< durat1 <<  "ms"  << std::endl;
        
        // 创建P1有限元实例
        femm::LagrangeP1FEM fem(mesh);

        std::__1::chrono::steady_clock::time_point now2 = std::chrono::steady_clock::now();
        double durat2 = std::chrono::duration_cast<std::chrono::milliseconds>(now2-now1).count();
        std::cout << "fem set! time: "<< durat2 <<  "ms" << std::endl;

        application::PoissonDirichlet poisson(fem);

        std::__1::chrono::steady_clock::time_point now3 = std::chrono::steady_clock::now();
        double durat3 = std::chrono::duration_cast<std::chrono::milliseconds>(now3-now2).count();
        std::cout << "PDE set! time: "<< durat3 <<  "ms" << std::endl;

        fempatch::FEMPatchDirichletZero dirichletBC(fem);
        
        
        PetscInt numnodes = fem.GetNumNodes();

        std::cout << "num of nodes: " << numnodes <<std::endl;

        poisson.AddBC(&dirichletBC);
        std::__1::chrono::steady_clock::time_point now4 = std::chrono::steady_clock::now();
        double durat4 = std::chrono::duration_cast<std::chrono::milliseconds>(now4-now3).count();
        std::cout << "BC applied! time: "<< durat4 <<  "ms" << std::endl;
        
        poisson.Prepare();
        std::__1::chrono::steady_clock::time_point now5 = std::chrono::steady_clock::now();
        double durat5 = std::chrono::duration_cast<std::chrono::milliseconds>(now5-now4).count();
        std::cout << "app prepared. time: "<<durat5 <<  "ms" << std::endl;
        
        poisson.Solve();
        std::__1::chrono::steady_clock::time_point now6 = std::chrono::steady_clock::now();
        double durat6 = std::chrono::duration_cast<std::chrono::milliseconds>(now6-now5).count();
        std::cout << "Solve Done! time: "<< durat6 <<  "ms" <<std::endl;

        double durattot = std::chrono::duration_cast<std::chrono::milliseconds>(now6-now1).count();
        std::cout << "total solve time: " <<durattot<< "ms" << std::endl;


        Vec &solution = poisson.GetSolverSol();
        Vec sol_exact;
        PetscCall(fem.DomainProject(constants::Exact,sol_exact));
        PetscScalar err;
        Mat stiff = poisson.GetSolverStiff();
        
        //PetscCall(numerical::VecErrL2Rel(solution,sol_exact,err));
        
        PetscCall(numerical::VecErrL2(solution,sol_exact,err));
        std::cout << "rel L2 err: " << err << std::endl;
        //PetscCall(VecView(sol_exact,PETSC_VIEWER_STDOUT_WORLD));
        //PetscCall(VecView(solution,PETSC_VIEWER_STDOUT_WORLD));

        return 0;
    }

}
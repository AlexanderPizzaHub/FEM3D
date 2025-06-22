#include "FEMapplication.hpp"

namespace application
{
    #pragma region poissondirichlet
    PoissonDirichlet::PoissonDirichlet(femm::LagrangeP1FEM &fem)
    {
        // 初始化FEM机器
        fem_ = &fem;
        std::cout << "PoissonDirichlet initialized!" << std::endl;
    };

    PoissonDirichlet::~PoissonDirichlet()
    {
        //delete fem_;
    };

    PetscErrorCode PoissonDirichlet::AddBC(fempatch::FEMPatchDirichlet *patch)
    {
        // 添加边界处理类
        patches_.push_back(patch);
        return 0;
    }

    PetscErrorCode PoissonDirichlet::Prepare()
    {
        // 组装矩阵和右端项
        if(patches_.empty())
        {
            std::cerr << "No boundary conditions added!" << std::endl;
            return -1; // 没有添加边界条件
        }

        PetscCall(fem_->AssembleStiff());

        Vec source_data;
        PetscCall(fem_->DomainProject(constants::Source,source_data));
        PetscCall(fem_->AssembleRHS(source_data));
        Mat stiff = fem_->GetStiff();
        Vec rhs = fem_->GetRHS();
        for (auto &patch : patches_)
        {
            PetscCall(patch->ApplyBC(stiff, stiff_sub, rhs, rhs_sub));
        }
        PetscCall(MatAssemblyBegin(stiff_sub, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff_sub, MAT_FINAL_ASSEMBLY));
        PetscCall(VecAssemblyBegin(rhs_sub));
        PetscCall(VecAssemblyEnd(rhs_sub));
        // 创建求解器
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
        PetscCall(KSPSetOperators(ksp, stiff_sub, stiff_sub));
        PetscCall(KSPSetFromOptions(ksp));

        // 创建解向量
        PetscCall(utils::VecSetup(fem_->GetNumNodes(), sol));
        PetscCall(VecDuplicate(rhs_sub, &sol_sub));
        return 0;
    }

    PetscErrorCode PoissonDirichlet::Solve()
    {
        // 解线性方程组
        
        PetscCall(KSPSolve(ksp, rhs_sub, sol_sub));
        for (auto &patch : patches_)
        {
            PetscInt patch_type = patch->GetPatchType();
            if(patch_type == 1)
            {
                PetscCall(patch->AugVec(sol_sub, sol));
            }
        }

        return 0;
    }
    #pragma endregion



    #pragma region poissonneumanndirichlet
    PoissonMixed::PoissonMixed(femm::LagrangeP1FEM &fem)
    {
        // 初始化FEM机器
        fem_ = &fem;
        std::cout << "PoissonMixed initialized!" << std::endl;
    };

    PoissonMixed::~PoissonMixed()
    {
        //delete fem_;
    };

    PetscErrorCode PoissonMixed::AddBCDirichlet(fempatch::FEMPatchDirichlet *patch)
    {
        // 添加边界处理类
        patches_dirichlet_.push_back(patch);
        return 0;
    }
    PetscErrorCode PoissonMixed::AddBCNeumann(fempatch::FEMPatchNeumann *patch)
    {
        // 添加边界处理类
        patches_neumann_.push_back(patch);
        return 0;
    }

    PetscErrorCode PoissonMixed::Prepare()
    {
        // 组装矩阵和右端项
        if(patches_dirichlet_.empty() & patches_neumann_.empty())
        {
            std::cerr << "No boundary conditions added!" << std::endl;
            return -1; // 没有添加边界条件
        }

        PetscCall(fem_->AssembleStiff());

        Vec source_data;
        PetscCall(fem_->DomainProject(constants::Source,source_data));
        PetscCall(fem_->AssembleRHS(source_data));
        Mat stiff = fem_->GetStiff();
        Vec rhs = fem_->GetRHS();


        for (auto &patch : patches_dirichlet_)
        {
            PetscCall(patch->ApplyBC(stiff, stiff_sub, rhs, rhs_sub));
        }

        for(auto &patch : patches_neumann_)
        {
            PetscCall(patch->ApplyBCBySource(source_data));
        }

        PetscCall(MatAssemblyBegin(stiff_sub, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(stiff_sub, MAT_FINAL_ASSEMBLY));
        PetscCall(VecAssemblyBegin(rhs_sub));
        PetscCall(VecAssemblyEnd(rhs_sub));
        // 创建求解器
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
        PetscCall(KSPSetOperators(ksp, stiff_sub, stiff_sub));
        PetscCall(KSPSetFromOptions(ksp));

        // 创建解向量
        PetscCall(utils::VecSetup(fem_->GetNumNodes(), sol));
        PetscCall(VecDuplicate(rhs_sub, &sol_sub));
        return 0;
    }

    PetscErrorCode PoissonMixed::Solve()
    {
        // 解线性方程组
        
        PetscCall(KSPSolve(ksp, rhs_sub, sol_sub));
        for (auto &patch : patches_dirichlet_)
        {
            PetscInt patch_type = patch->GetPatchType();
            if(patch_type == 1)
            {
                PetscCall(patch->AugVec(sol_sub, sol));
            }
        }

        return 0;
    }

    #pragma endregion
}
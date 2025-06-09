#include "FEMapplication.hpp"

namespace application
{
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

    PetscErrorCode PoissonDirichlet::AddBC(fempatch::FEMPatchBase *patch)
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
        PetscCall(fem_->AssembleRHS());
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
            PetscCall(patch->AugVec(sol_sub, sol));
        }

        return 0;
    }
}
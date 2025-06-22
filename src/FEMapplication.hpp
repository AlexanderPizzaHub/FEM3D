/*
此模块包含多个类，每个类负责完成一个特定形式PDE的求解。
*/
#pragma once
#include <petsc.h>
#include <string>
#include <vector>

#include "FEMmachine.hpp"
#include "FEMpatch.hpp"
#include "numericaltools.hpp"
#include "Const.hpp"


namespace application
{
    class PoissonDirichlet
    {
        public:
            // 暂时不考虑外部提供数据文件
            PoissonDirichlet(femm::LagrangeP1FEM &fem);
            ~PoissonDirichlet();

            PetscErrorCode AddBC(fempatch::FEMPatchDirichlet *patch);
        

            PetscErrorCode Prepare(); // 组矩阵，完成所有矩阵final assembly
            PetscErrorCode Solve(); // 解矩阵

            Mat &GetSolverStiff(){return stiff_sub;};
            Vec &GetSolverRHS(){return rhs_sub;};
            Vec &GetSolverSol(){return sol;};
            KSP &GetSolverKSP(){return ksp;};


        private:
            femm::LagrangeP1FEM *fem_; // 有限元方法数据结构
            std::vector<fempatch::FEMPatchDirichlet*> patches_; // 边界处理类

            // 如果用子矩阵法处理边界，那么必须有额外一个矩阵存储子矩阵，我们记为stiff_sub
            // 如果用罚方法，置一法等等不改变矩阵大小的方法，我们可以直接在原矩阵上做修改，以下私有变量可以存为原矩阵的引用
            // 我们暂时使用子矩阵法
            Mat stiff_sub;
            Vec rhs_sub;
            Vec sol_sub;
            Vec sol;
            KSP ksp;

    };


    class PoissonMixed
    {
        public:
            // 暂时不考虑外部提供数据文件
            PoissonMixed(femm::LagrangeP1FEM &fem);
            ~PoissonMixed();

            PetscErrorCode AddBCDirichlet(fempatch::FEMPatchDirichlet *patch);
            PetscErrorCode AddBCNeumann(fempatch::FEMPatchNeumann *patch);
        

            PetscErrorCode Prepare(); // 组矩阵，完成所有矩阵final assembly
            PetscErrorCode Solve(); // 解矩阵

            Mat &GetSolverStiff(){return stiff_sub;};
            Vec &GetSolverRHS(){return rhs_sub;};
            Vec &GetSolverSol(){return sol;};
            KSP &GetSolverKSP(){return ksp;};


        private:
            femm::LagrangeP1FEM *fem_; // 有限元方法数据结构
            std::vector<fempatch::FEMPatchDirichlet*> patches_dirichlet_; // 边界处理类
            std::vector<fempatch::FEMPatchNeumann*> patches_neumann_;

            // 如果用子矩阵法处理边界，那么必须有额外一个矩阵存储子矩阵，我们记为stiff_sub
            // 如果用罚方法，置一法等等不改变矩阵大小的方法，我们可以直接在原矩阵上做修改，以下私有变量可以存为原矩阵的引用
            // 我们暂时使用子矩阵法
            Mat stiff_sub;
            Vec rhs_sub;
            Vec sol_sub;
            Vec sol;
            KSP ksp;

    };
}
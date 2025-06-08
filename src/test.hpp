/*
此模块负责测试案例
*/
#pragma once 

#include <petsc.h>
#include <string>
#include <vector>
#include "mesh.hpp"
#include "FEMmachine.hpp"

namespace test
{
    PetscErrorCode TestMeshDMPlex();

    PetscErrorCode TestFEMMachine();
}
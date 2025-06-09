/*
此模块负责测试案例
*/
#pragma once 

#include <petsc.h>
#include <string>
#include <vector>
#include <iostream>
#include <chrono>

#include "mesh.hpp"
#include "FEMmachine.hpp"
#include "FEMpatch.hpp"
#include "FEMapplication.hpp"
namespace test
{
    PetscErrorCode TestMeshDMPlex();
    PetscErrorCode TestMeshVertMap();
    PetscErrorCode TestStiffDeg();

    PetscErrorCode TestFEMMachine();

    PetscErrorCode TestStiff();

    PetscErrorCode TestPatch();

    PetscErrorCode TestSolver();
}
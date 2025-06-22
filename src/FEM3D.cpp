/*
此模块放置main函数。
*/
char help[] = "FEM3D codes.";
#include <iostream>
#include <petsc.h>
#include "utils.hpp"
#include "fileIO.hpp"
#include "mesh.hpp"
#include "numericaltools.hpp"
#include "FEMmachine.hpp"
#include "FEMpatch.hpp"
#include "FEMapplication.hpp"
#include "test.hpp"

int main(int argc, char **argv)
{
    int ierr;
    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

     //test::TestMeshDMPlex();
     //test::TestMeshVertMap();
   test::TestFEMMachine();
    //test::TestStiff();
   // test::TestStiffDeg();
   //test::TestPatch();
   //test::TestSolver();
    
    PetscCall(PetscFinalize());
    return 0;
}
#include "mesh.hpp"
#include <iostream>

namespace mesh
{
    MeshDMPlex::MeshDMPlex(const char filename[]){
        DMPlexCreateGmshFromFile(PETSC_COMM_WORLD, filename, PETSC_TRUE, &dm_);
        MarkCubeBoundary(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    };

    // 第三种读取方式，我们手动给DMPlex标记边界。网格边界应该由mesh文件指定，因此该功能计划在未来废弃。
    //MeshDMPlex::MeshDMPlex(const char filename[], const std::string domain_type);

    MeshDMPlex::~MeshDMPlex(){};

    // 用于手动标记边界的函数。计划在未来废弃。
    PetscErrorCode MeshDMPlex::MarkCubeBoundary( PetscScalar left, PetscScalar right, PetscScalar front, PetscScalar back, PetscScalar top, PetscScalar bottom){
        using namespace numerical;
        // 假设msh中各维度几何实体的序号连续排列
        PetscInt node_start, node_end, p;
        PetscScalar *global_coords, *coords;
        Vec coordinates;

        PetscCall(DMLabelCreate(PETSC_COMM_WORLD, "boundary", &label_));

        PetscCall(DMPlexGetDepthStratum(dm_, 0, &node_start, &node_end));
        PetscCall(DMGetCoordinatesLocal(dm_, &coordinates));

        VecGetArray(coordinates, &global_coords);

        coords = new PetscScalar[dim];
        for (p = node_start; p < node_end; ++p)
        {
            //std::cout << "Processing point: " << p-node_start << std::endl;
            coords[0] = global_coords[dim * (p - node_start)];
            coords[1] = global_coords[dim * (p - node_start) + 1];
            coords[2] = global_coords[dim * (p - node_start) + 2];
            //std::cout << "Coordinates: (" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
            if (
                near(coords[0], left) ||
                near(coords[0], right) || 
                near(coords[1], front) || 
                near(coords[1], back) || 
                near(coords[2], top) || 
                near(coords[2], bottom))
            {
                //std::cout << "Marking boundary point: " << p << " at coords: (" 
                          //<< coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
                PetscCall(DMLabelSetValue(label_, p, 1));
            }
            else
            {
                //std::cout << "Marking inner point: " << p << " at coords: (" 
                          //<< coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
                PetscCall(DMLabelSetValue(label_, p, 0));
            }
        }
        PetscCall(DMPlexLabelComplete(dm_, label_));

        //PetscCall(DMLabelGetStratumIS(label_, 1, &bdryIS_));
        //PetscCall(DMLabelGetStratumIS(label_, 0, &innerIS_));
        PetscCall(GetISbyLabel(1, bdryIS_));
        PetscCall(GetISbyLabel(0, innerIS_));
        PetscCall(ISShift(innerIS_, -node_start, innerIS_));
        PetscCall(ISShift(bdryIS_, -node_start, bdryIS_));

        return 0;
    };


    // 根据label获取对应的网格点序号。
    PetscErrorCode MeshDMPlex::GetISbyLabel(PetscInt label_index, IS &is)
    {
        PetscCall(DMLabelGetStratumIS(label_, label_index, &is));
        return 0;
    };

}
#include "mesh.hpp"
#include <iostream>

namespace mesh
{
    MeshDMPlex::MeshDMPlex(const char filename[]){
        DMPlexCreateGmshFromFile(PETSC_COMM_WORLD, filename, PETSC_TRUE, &dm_);
       //MarkCubeBoundary(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
       MarkMixBoundary(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        MakeCell2VertMap();
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

    PetscErrorCode MeshDMPlex::MarkMixBoundary( PetscScalar left, PetscScalar right, PetscScalar front, PetscScalar back, PetscScalar top, PetscScalar bottom){
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
                near(coords[1], front) || 
                near(coords[1], back) 
                )
            {
                //std::cout << "Marking boundary point: " << p << " at coords: (" 
                          //<< coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
                PetscCall(DMLabelSetValue(label_, p, 1));
            }
            else if
            (
                near(coords[0], left) ||
                near(coords[0], right) || 
                near(coords[2], top) || 
                near(coords[2], bottom)
            )
            {
                PetscCall(DMLabelSetValue(label_, p, 2));
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
        //PetscCall(GetISbyLabel(1, bdryIS_));
        //PetscCall(GetISbyLabel(0, innerIS_));
        //PetscCall(ISShift(innerIS_, -node_start, innerIS_));
        //PetscCall(ISShift(bdryIS_, -node_start, bdryIS_));

        return 0;
    };


    // 根据label获取对应的网格点序号。
    PetscErrorCode MeshDMPlex::GetISbyLabel(PetscInt label_index, IS &is)
    {
        PetscCall(DMLabelGetStratumIS(label_, label_index, &is));
        return 0;
    };

    PetscErrorCode MeshDMPlex::GetFace2VertIdx(PetscInt face_idx, PetscInt* face2node_idx)
    {
        const PetscInt *edgecone, *nodecone;
        PetscInt node_start, node_end, edge_conesize, node_conesize;
        PetscInt neighbor_count;
        PetscInt *node_neighbors_tmp;

        PetscCall(DMPlexGetDepthStratum(dm_, 0, &node_start, &node_end));
        PetscCall(DMPlexGetCone(dm_, face_idx, &edgecone));
        PetscCall(DMPlexGetConeSize(dm_, face_idx, &edge_conesize));

        assert(edge_conesize == 3);

        node_neighbors_tmp = new PetscInt[6];
        neighbor_count = 0;
        for (PetscInt i = 0; i < edge_conesize; i++)
        {
            PetscInt edge = edgecone[i];
            PetscCall(DMPlexGetConeSize(dm_, edge, &node_conesize));
            PetscCall(DMPlexGetCone(dm_, edge, &nodecone));

            assert(node_conesize == 2);
            for (PetscInt j = 0; j < node_conesize; j++)
            {
                PetscInt node = nodecone[j];
                assert(node >= node_start && node < node_end);
                node_neighbors_tmp[neighbor_count] = node;
                neighbor_count++;
            }
        }
        assert(neighbor_count == 6);
        PetscCall(PetscSortInt(neighbor_count, node_neighbors_tmp));

        for (PetscInt i = 0; i < 3; i++)
        {
            face2node_idx[i] = node_neighbors_tmp[2 * i];
        }
        return 0;
    }

    PetscErrorCode MeshDMPlex::GetCell2VertIdx(PetscInt cell_idx, PetscInt* cell2node_idx)
    {
        // 重点检查
        const PetscInt *facecone, *nodecone;
        PetscInt node_start, node_end, face_conesize, node_conesize;
        PetscInt neighbor_count;
        PetscInt *node_neighbors_tmp, *node_neighbors_sub;

        PetscCall(DMPlexGetDepthStratum(dm_, 0, &node_start, &node_end));
        PetscCall(DMPlexGetCone(dm_, cell_idx, &facecone));
        PetscCall(DMPlexGetConeSize(dm_, cell_idx, &face_conesize));
        //std::cout << "cell_idx: " << cell_idx << ", face_conesize: " << face_conesize << std::endl;
        assert(face_conesize == 4);
        

        node_neighbors_tmp = new PetscInt[12]; // 4面*3点
        node_neighbors_sub = new PetscInt[3]; // 3个点
        neighbor_count = 0;
       
        for (PetscInt i = 0; i < face_conesize; i++)
        {
            PetscInt face = facecone[i];
            PetscCall(GetFace2VertIdx(face, node_neighbors_sub));
            //std::cout << "here!!" << std::endl;

            for(PetscInt j = 0; j < 3; j++)
            {
                node_neighbors_tmp[3*i + j] = node_neighbors_sub[j];
                neighbor_count++;
            }
        }
        assert(neighbor_count == 12);
        

        PetscCall(PetscSortInt(neighbor_count, node_neighbors_tmp));

        for (PetscInt i = 0; i < 4; i++)
        {
            cell2node_idx[i] = node_neighbors_tmp[3 * i];
            //std::cout << "cell2node_idx[" << i << "]: " << cell2node_idx[i] << std::endl;
        }
        
        return 0;
    };

    PetscErrorCode MeshDMPlex::MakeCell2VertMap()
    {
        PetscInt cell_start, cell_end;
        PetscInt *cell2vert_idx;

        // !!!!!!!!!!!!!
        PetscCall(DMPlexGetDepthStratum(dm_, dim, &cell_start, &cell_end));
        cell2vert_map_.reserve(cell_end-cell_start);
        cell2vert_idx = new PetscInt[4]; // 四面体的四个顶点

        for (PetscInt p_cell = cell_start; p_cell < cell_end; ++p_cell)
        {
            PetscCall(GetCell2VertIdx(p_cell,cell2vert_idx));
            std::vector<PetscInt> vert_indices(4);
            for (PetscInt i = 0; i < 4; i++)
            {
                vert_indices[i] = cell2vert_idx[i];
            }
            cell2vert_map_.push_back(vert_indices);
        }
        return 0;
    };

}
//
// Created by Hallison da Paz on 15/05/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H
#define DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

#include <vector>
#include <list>
#include "glm/glm.hpp"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "old/qef.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  DefaultMesh;

struct Vertex {
    // position coordinates
    glm::vec3 position;
    // normal coordinates
    glm::vec3 normal;
    // texture coordinates
    //glm::vec2 texCoords;
    //color coordinates
    glm::uvec3 color;

    Vertex(glm::vec3 pos):
            position(pos),
            normal(0),
            //texCoords(0)
            color(0){
    }

    Vertex(): position(0), normal(0), color(0){}
};

struct Triangle{
    // Triangle vertices
    int a;
    int b;
    int c;
};

typedef std::vector<Vertex> VertexBuffer;
typedef std::vector<Triangle> IndexBuffer;

struct OctreeMeshInfo {
// ---------------------------------------------------------------------------- OPENMESH
    std::list<DefaultMesh::FaceHandle> innerFaces;
    std::list<DefaultMesh::FaceHandle> crossingFaces;
};

struct OctreeDrawInfo
{
    OctreeDrawInfo()
            : index(-1)
            , corners(0)
    {
    }

    int	index;
    int	corners;
    glm::vec3 position;
    glm::vec3 averageNormal;
    svd::QefData qef;
};

#endif //DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

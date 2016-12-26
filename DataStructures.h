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
#include "qef.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  DefaultMesh;

//glm::uvec3 defaultgray = glm::uvec3(128, 128, 128); //middle gray

struct Vertex {
    // position coordinates
    glm::vec3 position;
    // normal coordinates
    glm::vec3 normal;
    //color coordinates
    glm::uvec3 color;

    Vertex(glm::vec3 pos): position(pos), normal(0), color(glm::uvec3(128, 128, 128)){}

    Vertex(): position(0), normal(0), color(glm::uvec3(128, 128, 128)){}
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

// ----------------------------------------------------------------------------

enum OctreeNodeType
{
    NODE_NONE,
    NODE_INTERNAL,
    NODE_PSEUDO,
    NODE_LEAF,
};

// ----------------------------------------------------------------------------

#endif //DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

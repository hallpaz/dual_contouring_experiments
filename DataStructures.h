//
// Created by Hallison da Paz on 15/05/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H
#define DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

#include <vector>
#include "glm/glm.hpp"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <list>
#include "Constants.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  DefaultMesh;

struct Vertex {
    // position coordinates
    vecr position;
    // normal coordinates
    vecr normal;
    // texture coordinates
    //glm::vec2 texCoords;
    //color coordinates
    glm::uvec3 color;

    Vertex(vecr pos):
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

struct OctreeMeshInfo {
// ---------------------------------------------------------------------------- OPENMESH
    std::list<DefaultMesh::FaceHandle> innerFaces;
    std::list<DefaultMesh::FaceHandle> crossingFaces;
};

typedef std::vector<Vertex> VertexBuffer;
typedef std::vector<Triangle> IndexBuffer;


#endif //DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

//
// Created by Hallison da Paz on 15/05/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H
#define DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

#include <vector>
#include "glm/glm.hpp"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<>  DefaultMesh;

struct Vertex {
    // position coordinates
    glm::vec3 position;
    // normal coordinates
    glm::vec3 normal;
    // texture coordinates
    glm::vec2 texCoords;

    Vertex(glm::vec3 pos):
            position(pos),
            normal(0),
            texCoords(0){
    }
};

struct Triangle{
    // Triangle vertices
    int a;
    int b;
    int c;
};

typedef std::vector<Vertex> VertexBuffer;
typedef std::vector<Triangle> IndexBuffer;


#endif //DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

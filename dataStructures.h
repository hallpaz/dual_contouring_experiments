//
// Created by Hallison da Paz on 15/05/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H
#define DUAL_CONTOURING_EXPERIMENTS_DATASTRUCTURES_H

#include <vector>
#include "glm/glm.hpp"


const unsigned int posOffset3D = 0;
const unsigned int normalOffset3D = 3*sizeof(float);
const unsigned int texOffset3D= 6*sizeof(float);

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

struct Vertex2D {
    // position coordinates
    glm::vec2 position;
    // normal coordinates
    glm::vec3 normal;
    // texture coordinates
    glm::vec2 texCoords;
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

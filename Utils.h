//
//  Utils.hpp
//  PanoramicTour
//
//  Created by Hallison da Paz on 06/11/2015.
//  Copyright © 2015 Visgraf. All rights reserved.
//

#ifndef Utils_hpp
#define Utils_hpp

#include <vector>
#include <string>

#include "glm/glm.hpp"
#include "dataStructures.h"


void write_Ply(std::vector<Vertex> &vertices, std::vector<Triangle> &faces, std::string filename);
void write_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces);

float read_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces, glm::vec3 &min);

std::string hashvertex(const glm::vec3& vertex);
std::string hashedge(const glm::vec3& v1, const glm::vec3& v2);

void void barycentric(glm::vec3 p, glm::vec3 a, glm::vec3 b, glm::vec3 c, float &u, float &v, float &w);

#endif /* Utils_hpp */

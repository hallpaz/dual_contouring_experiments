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
void write_OFF(std::vector<Vertex> &vertices, std::vector<Triangle> &faces, std::string filename);

float read_OFF(std::vector<Vertex> &vertices, std::vector<Triangle> &faces, std::string filename);

std::string hashvertex(const glm::vec3& vertex);
std::string hashedge(const glm::vec3& v1, const glm::vec3& v2);

#endif /* Utils_hpp */

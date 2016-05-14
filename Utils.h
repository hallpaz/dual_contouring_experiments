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

void write_Ply(std::vector<glm::vec3> &vertices, std::vector<int> &faces, std::string filename);
void write_OFF(std::vector<glm::vec3> &vertices, std::vector<int> &faces, std::string filename);

#endif /* Utils_hpp */

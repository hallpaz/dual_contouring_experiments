//
// Created by Hallison da Paz on 03/10/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_RECONSTRUCTION_H
#define DUAL_CONTOURING_EXPERIMENTS_RECONSTRUCTION_H

#include "octree.h"

namespace Fusion
{
    OctreeNode* octree_from_samples(const glm::vec3 &min, const float size, const int height, std::vector<std::string> meshfiles);
}


#endif //DUAL_CONTOURING_EXPERIMENTS_RECONSTRUCTION_H

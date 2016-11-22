//
// Created by Hallison da Paz on 03/10/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_RECONSTRUCTION_H
#define DUAL_CONTOURING_EXPERIMENTS_RECONSTRUCTION_H

#include "Octree.h"

namespace Fusion
{
    OctreeNode *octree_from_samples(const glm::vec3 &min, const float size, const unsigned int max_depth,
                                        std::vector<std::string> meshfiles, std::vector<glm::vec3> cameras);
}


#endif //DUAL_CONTOURING_EXPERIMENTS_RECONSTRUCTION_H

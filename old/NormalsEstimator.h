//
// Created by hallpaz on 9/1/16.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_NORMALSESTIMATOR_H
#define DUAL_CONTOURING_EXPERIMENTS_NORMALSESTIMATOR_H

#include "../DataStructures.h"
// -------------------------------------------

#include <string>

enum Quality {
    COARSE, FINE
};

class NormalsEstimator
{
public:
    static void update_model_normals(std::string model_path, std::string output_path, Quality q = FINE);
    static void compute_better_normals(DefaultMesh &mesh);

};


#endif //DUAL_CONTOURING_EXPERIMENTS_NORMALSESTIMATOR_H

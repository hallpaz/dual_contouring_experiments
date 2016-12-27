//
// Created by Hallison da Paz on 27/12/2016.
//

#ifndef DUAL_CONTOURING_SANDBOX_TAOJUQEF_H
#define DUAL_CONTOURING_SANDBOX_TAOJUQEF_H


#include "GeoCommon.h"

#define TOLERANCE 0.0001f

class TaojuQef {
public:
    static float calcPoint ( float halfA[], float b[], float btb, float midpoint[], float rvalue[], BoundingBoxf &box, float *mat );
};


#endif //DUAL_CONTOURING_SANDBOX_TAOJUQEF_H

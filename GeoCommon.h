//
// Created by Hallison da Paz on 27/12/2016.
//
/*
  Structure definitions for points and triangles.
  Copyright (C) 2011  Tao Ju
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef DUAL_CONTOURING_SANDBOX_GEOCOMMON_H
#define DUAL_CONTOURING_SANDBOX_GEOCOMMON_H


#define USE_MINIMIZER

#include "glm/vec3.hpp"


typedef struct {
    glm::vec3 begin;
    glm::vec3 end;
} BoundingBoxf;

#endif //DUAL_CONTOURING_SANDBOX_GEOCOMMON_H

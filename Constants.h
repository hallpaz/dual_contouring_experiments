//
// File Created by Hallison da Paz on 02/09/2016.
// Constants extracted from Tao Ju implementation
/*
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


#ifndef DUAL_CONTOURING_EXPERIMENTS_CONSTANTS_H
#define DUAL_CONTOURING_EXPERIMENTS_CONSTANTS_H

#include "glm/glm.hpp"


typedef float Real;
// ----------------------------------------------------------------------------

const int MATERIAL_AIR = 1;
const int MATERIAL_SOLID = 0;
const int MATERIAL_UNKNOWN = -1;
const int MATERIAL_AMBIGUOUS = -2;
const unsigned short NUM_CHILDREN = 8;
const unsigned short NUM_EDGES = 12;
const Real EPSILON = 1e-5f;

const float QEF_ERROR = 1e-2f;
const int QEF_SWEEPS = 4;
// ----------------------------------------------------------------------------
const int vneighbors[8][3] = {
        {1, 2, 4},
        {0, 3, 5},
        {0, 3, 6},
        {1, 2, 7},
        {0, 5, 6},
        {1, 4, 7},
        {2, 4, 7},
        {3, 5, 6}
};
// ----------------------------------------------------------------------------

const glm::vec3 CHILD_MIN_OFFSETS[] =
        {
                // needs to match the vertMap from Dual Contouring impl
                glm::vec3( 0, 0, 0 ),
                glm::vec3( 0, 0, 1 ),
                glm::vec3( 0, 1, 0 ),
                glm::vec3( 0, 1, 1 ),
                glm::vec3( 1, 0, 0 ),
                glm::vec3( 1, 0, 1 ),
                glm::vec3( 1, 1, 0 ),
                glm::vec3( 1, 1, 1 ),
        };

// ----------------------------------------------------------------------------
// data from the original DC impl, drives the contouring process

const int edgevmap[12][2] =
        {
                {0,4},{1,5},{2,6},{3,7},	// x-axis
                {0,2},{1,3},{4,6},{5,7},	// y-axis
                {0,1},{2,3},{4,5},{6,7}		// z-axis
        };

const int edgemask[3] = { 5, 3, 6 } ;

const int vertMap[8][3] =
        {
                {0,0,0},
                {0,0,1},
                {0,1,0},
                {0,1,1},
                {1,0,0},
                {1,0,1},
                {1,1,0},
                {1,1,1}
        };

const int faceMap[6][4] = {{4, 8, 5, 9}, {6, 10, 7, 11},{0, 8, 1, 10},{2, 9, 3, 11},{0, 4, 2, 6},{1, 5, 3, 7}} ;
const int cellProcFaceMask[12][3] = {{0,4,0},{1,5,0},{2,6,0},{3,7,0},{0,2,1},{4,6,1},{1,3,1},{5,7,1},{0,1,2},{2,3,2},{4,5,2},{6,7,2}} ;
const int cellProcEdgeMask[6][5] = {{0,1,2,3,0},{4,5,6,7,0},{0,4,1,5,1},{2,6,3,7,1},{0,2,4,6,2},{1,3,5,7,2}} ;

const int faceProcFaceMask[3][4][3] = {
        {{4,0,0},{5,1,0},{6,2,0},{7,3,0}},
        {{2,0,1},{6,4,1},{3,1,1},{7,5,1}},
        {{1,0,2},{3,2,2},{5,4,2},{7,6,2}}
} ;

const int faceProcEdgeMask[3][4][6] = {
        {{1,4,0,5,1,1},{1,6,2,7,3,1},{0,4,6,0,2,2},{0,5,7,1,3,2}},
        {{0,2,3,0,1,0},{0,6,7,4,5,0},{1,2,0,6,4,2},{1,3,1,7,5,2}},
        {{1,1,0,3,2,0},{1,5,4,7,6,0},{0,1,5,0,4,1},{0,3,7,2,6,1}}
};

const int edgeProcEdgeMask[3][2][5] = {
        {{3,2,1,0,0},{7,6,5,4,0}},
        {{5,1,4,0,1},{7,3,6,2,1}},
        {{6,4,2,0,2},{7,5,3,1,2}},
};

const int processEdgeMask[3][4] = {{3,2,1,0},{7,5,6,4},{11,10,9,8}} ;

// -------------------------------------------------------------------------------


const float POINT_DISTANCE_THRESHOLD = 0.000001f;
#endif //DUAL_CONTOURING_EXPERIMENTS_CONSTANTS_H

//if (mesh.is_boundary(*face))
//{
//for (auto fh_iter = mesh.cfh_iter(*face); fh_iter.is_valid(); ++fh_iter)
//{
//DefaultMesh::Normal faceNormal = mesh.normal(*face);
//DefaultMesh::Point middle_point = (mesh.point(mesh.to_vertex_handle(*fh_iter)) + mesh.point(mesh.from_vertex_handle((*fh_iter))))/2;;
//DefaultMesh::Point edge = mesh.point(mesh.to_vertex_handle(*fh_iter)) - mesh.point(mesh.from_vertex_handle((*fh_iter)));
//DefaultMesh::Normal planeNormal;
//if (mesh.is_boundary(*fh_iter))
//{
//planeNormal =  edge % faceNormal;
//planeNormal.normalize();
//featureQef.add(middle_point[0], middle_point[1], middle_point[2], planeNormal[0], planeNormal[1], planeNormal[2]);
//}
//DefaultMesh::HalfedgeHandle opposite_halfedge = mesh.opposite_halfedge_handle(*fh_iter);
//if (mesh.is_boundary(opposite_halfedge))
//{
//planeNormal = faceNormal % edge;
//planeNormal.normalize();
//featureQef.add(middle_point[0], middle_point[1], middle_point[2], planeNormal[0], planeNormal[1], planeNormal[2]);
//}
//}
//}
//                for (auto fh_iter = mesh.cfh_iter(*face); fh_iter.is_valid(); ++fh_iter)
//                {
//                    DefaultMesh::Normal faceNormal;
//                    DefaultMesh::Point middle_point = (mesh.point(mesh.to_vertex_handle(*fh_iter)) + mesh.point(mesh.from_vertex_handle((*fh_iter))))/2;;
//                    DefaultMesh::Point edge = mesh.point(mesh.to_vertex_handle(*fh_iter)) - mesh.point(mesh.from_vertex_handle((*fh_iter)));
//                    DefaultMesh::Normal planeNormal;
//                    if (!mesh.is_boundary(*fh_iter))
//                    {
//                        faceNormal = mesh.normal(*face);
//                        planeNormal =  edge % faceNormal;
//                        planeNormal.normalize();
//                        featureQef.add(middle_point[0], middle_point[1], middle_point[2], planeNormal[0], planeNormal[1], planeNormal[2]);
//                    }
//                    DefaultMesh::HalfedgeHandle opposite_halfedge = mesh.opposite_halfedge_handle(*fh_iter);
//                    if (!mesh.is_boundary(opposite_halfedge))
//                    {
//                        faceNormal = mesh.normal(mesh.face_handle(opposite_halfedge));
//                        planeNormal = faceNormal % edge;
//                        planeNormal.normalize();
//                        featureQef.add(middle_point[0], middle_point[1], middle_point[2], planeNormal[0], planeNormal[1], planeNormal[2]);
//                    }
//                }
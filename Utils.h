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
#include <unordered_map>
#include <list>

#include "glm/glm.hpp"
#include "DataStructures.h"

#include "Constants.h"
// ----------------------------------------------------------------------------
enum RelativePosition
{
    INSIDE,
    CROSSING,
    OUTSIDE,
};

class OctreeNode;
// ----------------------------------------------------------------------------

void write_Ply(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces);
void write_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces);
Real read_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces, vecr &min);

std::string hashvertex(const vecr &vertex);
std::string hashedge(const vecr &v1, const vecr &v2);

void barycentric(vecr p, vecr a, vecr b, vecr c, Real &u, Real &v, Real &w);
bool moller_triangle_intersection(vecr v1, vecr v2, Vertex *triangle_vertices, vecr &intersection_point);

vecr ApproximateZeroCrossingPosition(const vecr &p0, const vecr &p1);
vecr CalculateSurfaceNormal(const vecr &p);

RelativePosition vertexRelativePosition(const DefaultMesh &mesh, const DefaultMesh::VertexHandle &vertexHandle, vecr min, Real size);
RelativePosition halfedgeRelativePosition(const DefaultMesh &mesh, const DefaultMesh::HalfedgeHandle &halfedgeHandle, vecr min, Real size);
RelativePosition triangleRelativePosition(const DefaultMesh &mesh, const DefaultMesh::FaceHandle &faceHandle, vecr min, Real size);
RelativePosition vertexRelativePosition(const Vertex &vertex, vecr min, Real size);
RelativePosition triangleRelativePosition(const Vertex &a, const Vertex& b, const Vertex& c, vecr min, Real size);
RelativePosition vertexRelativePosition(const Vertex &vertex, vecr min, Real size);
RelativePosition triangleRelativePosition(const Vertex &a, const Vertex &b, const Vertex &c, vecr min, Real size);

vecr openmesh_to_glm(const OpenMesh::VectorT<float, 3> om_vec);


int computeSideOfPoint(const vecr point, const vecr intersection, const vecr face_normal);

int ray_mesh_intersection(const vecr origin, const vecr dest, DefaultMesh &mesh);

void updateVertexpool(std::unordered_map<std::string, int> &pool, const vecr &vertex, int &sign);

void updateSignsArray(int *vecsigns, int size);

void divideFacesByLocation(OctreeNode *node, std::list<DefaultMesh::FaceHandle> &facesList, const DefaultMesh &mesh);

// ----------------------------------------------------------------------------
inline void trace(std::string debug_text)
{
    std::cout << debug_text << std::endl;
}

// ----------------------------------------------------------------------------
void select_inner_crossing_faces(OctreeNode *node, const DefaultMesh &mesh);

#endif /* Utils_hpp */
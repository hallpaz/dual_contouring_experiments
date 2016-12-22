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


// ----------------------------------------------------------------------------
enum RelativePosition
{
    INSIDE,
    CROSSING,
    OUTSIDE,
};

class OctreeNode;
// ----------------------------------------------------------------------------

void write_Ply(std::string filename, VertexBuffer &vertices, IndexBuffer &faces);
void write_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces);
float read_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces, glm::vec3 &min);

std::string hashvertex(const glm::vec3 &vertex);
std::string hashedge(const glm::vec3 &v1, const glm::vec3 &v2);

void barycentric(glm::vec3 p, glm::vec3 a, glm::vec3 b, glm::vec3 c, float &u, float &v, float &w);
bool moller_triangle_intersection(glm::vec3 v1, glm::vec3 v2, Vertex *triangle_vertices, glm::vec3 &intersection_point);

glm::vec3 ApproximateZeroCrossingPosition(const glm::vec3 &p0, const glm::vec3 &p1);
glm::vec3 CalculateSurfaceNormal(const glm::vec3 &p);

RelativePosition vertexRelativePosition(const DefaultMesh &mesh, const DefaultMesh::VertexHandle &vertexHandle, glm::vec3 min, float size);
RelativePosition halfedgeRelativePosition(const DefaultMesh &mesh, const DefaultMesh::HalfedgeHandle &halfedgeHandle, glm::vec3 min, float size);
RelativePosition triangleRelativePosition(const DefaultMesh &mesh, const DefaultMesh::FaceHandle &faceHandle, glm::vec3 min, float size);
RelativePosition vertexRelativePosition(const Vertex &vertex, glm::vec3 min, float size);
RelativePosition triangleRelativePosition(const Vertex &a, const Vertex& b, const Vertex& c, glm::vec3 min, float size);
RelativePosition vertexRelativePosition(const Vertex &vertex, glm::vec3 min, float size);
RelativePosition triangleRelativePosition(const Vertex &a, const Vertex &b, const Vertex &c, glm::vec3 min, float size);

glm::vec3 openmesh_to_glm(const OpenMesh::VectorT<float, 3> om_vec);


int computeSideOfPoint(const glm::vec3 point, const glm::vec3 intersection, const glm::vec3 face_normal);

int ray_mesh_intersection(const glm::vec3 origin, const glm::vec3 dest, DefaultMesh &mesh);

void updateVertexpool(std::unordered_map<std::string, int> &pool, const glm::vec3 &vertex, int &sign);

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
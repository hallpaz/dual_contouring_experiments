/*

Implementations of Octree member functions.

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

#include <fstream>

#include "Constants.h"
#include "Contouring.h"
#include "density.h"
#include "octree.h"
#include "Utils.h"

#include "glm/ext.hpp"


using glm::vec3;
using glm::ivec3;


// -------------------------------------------------------------------------------

std::unordered_map<std::string, int> OctreeNode::vertexpool;
std::unordered_map<std::string, HermiteData> OctreeNode::edgepool;

// ----------------------------------------------------------------------------

int computeSideOfPoint(const glm::vec3 point, const glm::vec3 intersection, const glm::vec3 face_normal)
{
    return glm::dot(point - intersection, face_normal) < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
}

// ----------------------------------------------------------------------------
void updateVertexpool(std::unordered_map<std::string, int> &pool, const glm::vec3 &vertex, int &sign)
{
    if(sign == MATERIAL_UNKNOWN)
    {
        std::cout << "TENTA OUTRA NEGAO!!" << std::endl;
        exit(98);
    }
    std::string vertexhash = hashvertex(vertex);
    if (pool.count(vertexhash) == 0)
    {
        pool[vertexhash] = sign;
    }
    else
    {
        if (pool[vertexhash] == MATERIAL_SOLID)
        {   // if it was marked as interior anytime
            sign = pool[vertexhash];
        }
        else
        {
                pool[vertexhash] = sign;
        }
    }
    if (sign != MATERIAL_UNKNOWN && pool[vertexhash] == MATERIAL_UNKNOWN){
        std::cout << "LEFTING UNKNOWN!!!!!!!" << std::endl;
        exit(98);
    }
}

// ----------------------------------------------------------------------------
void updateSignsArray(int *vecsigns, int size)
{
    bool checksigns = true;
    while(checksigns) {
        checksigns = false;
        for (size_t i = 0; i < size; ++i) {
            if (vecsigns[i] == MATERIAL_UNKNOWN) {
                checksigns = true;
                for (int j = 0; j < 3; ++j) {
                    int n = vneighbors[i][j];
                    if (vecsigns[n] != MATERIAL_UNKNOWN) {
                        vecsigns[i] = vecsigns[n];
                        break;
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------

OctreeNode* ConstructLeaf(OctreeNode* leaf)
{

    if (!leaf || leaf->height != 0)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }

    int corners = 0;
    for (int i = 0; i < 8; i++)
    {
        const vec3 cornerPos = leaf->min + CHILD_MIN_OFFSETS[i]*leaf->size;
        const float density = Density_Func(vec3(cornerPos));
        const int material = density < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
        corners |= (material << i);
    }

    if (corners == 0 || corners == 255)
    {
        // voxel is full inside or outside the volume
        delete leaf;
        return nullptr;
    }

    // otherwise the voxel contains the surface, so find the edge intersections
    const int MAX_CROSSINGS = 6;
    int edgeCount = 0;
    vec3 averageNormal(0.f);
    svd::QefSolver qef;

    for (int i = 0; i < 12 && edgeCount < MAX_CROSSINGS; i++)
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];

        //take the signal of each vertex of the edge
        const int m1 = (corners >> c1) & 1;
        const int m2 = (corners >> c2) & 1;

        if ((m1 == MATERIAL_AIR && m2 == MATERIAL_AIR) ||
            (m1 == MATERIAL_SOLID && m2 == MATERIAL_SOLID))
        {
            // no zero crossing on this edge
            continue;
        }

        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);
        const vec3 p = ApproximateZeroCrossingPosition(p1, p2);
        const vec3 n = CalculateSurfaceNormal(p);
        qef.add(p.x, p.y, p.z, n.x, n.y, n.z);

        averageNormal += n;

        edgeCount++;
    }

    svd::Vec3 qefPosition;
    qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);

    OctreeDrawInfo* drawInfo = new OctreeDrawInfo;
    drawInfo->position = vec3(qefPosition.x, qefPosition.y, qefPosition.z);
    drawInfo->qef = qef.getData();

    const vec3 min = vec3(leaf->min);
    const vec3 max = vec3(leaf->min + vec3(leaf->size));
    if (drawInfo->position.x < min.x || drawInfo->position.x > max.x ||
        drawInfo->position.y < min.y || drawInfo->position.y > max.y ||
        drawInfo->position.z < min.z || drawInfo->position.z > max.z)
    {
        const auto& mp = qef.getMassPoint();
        drawInfo->position = vec3(mp.x, mp.y, mp.z);
    }

    drawInfo->averageNormal = glm::normalize(averageNormal / (float)edgeCount);
    drawInfo->corners = corners;

    leaf->type = Node_Leaf;
    leaf->drawInfo = drawInfo;

    return leaf;
}

// -------------------------------------------------------------------------------

OctreeNode* ConstructOctreeNodes(OctreeNode* node)
{
    if (!node)
    {
        std::cout << "Trying to construct empty node" << std::endl;
        return nullptr;
    }

    if (node->height == 0)
    {
        return ConstructLeaf(node);
    }

    const float childSize = node->size / 2;
    const int childHeight = node->height - 1;
    bool hasChildren = false;

    for (int i = 0; i < 8; i++)
    {
        OctreeNode* child = new OctreeNode;
        child->size = childSize;
        child->height = childHeight;
        child->min = node->min + (CHILD_MIN_OFFSETS[i] * childSize);
        child->type = Node_Internal;

        node->children[i] = ConstructOctreeNodes(child);
        hasChildren |= (node->children[i] != nullptr);
    }

    if (!hasChildren)
    {
        delete node;
        return nullptr;
    }

    return node;
}

// -------------------------------------------------------------------------------

OctreeNode* BuildOctree(const vec3& min, const float size, const int height, const float threshold)
{
    OctreeNode* root = new OctreeNode;
    root->min = min;
    root->size = size;
    root->height = height;
    root->type = Node_Internal;

    std::cout << "Octree.cpp: will construct nodes" << std::endl;
    ConstructOctreeNodes(root);
    std::cout << "Octree.cpp: will simplify nodes" << std::endl;
    root = SimplifyOctree(root, threshold);
    std::cout << "Octree.cpp: did simplify nodes" << std::endl;

    return root;
}

// -------------------------------------------------------------------------------

void DestroyOctree(OctreeNode* node)
{
    if (!node)
    {
        return;
    }

    for (int i = 0; i < 8; i++)
    {
        DestroyOctree(node->children[i]);
    }

    if (node->drawInfo)
    {
        delete node->drawInfo;
    }

    delete node;
}

// -------------------------------------------------------------------------------

OctreeNode* BuildOctreeFromOpenMesh(const glm::vec3 &min, const float size, const int height, const DefaultMesh &mesh)
{
    OctreeNode* root = new OctreeNode;
    root->min = min;
    root->size = size;
    root->height = height;
    root->type = Node_Internal;

    std::cout << "Octree.cpp: will construct nodes" << std::endl;
    ConstructOctreeNodesFromOpenMesh(root, mesh);
    std::cout << "Octree is Ready" << std::endl;
    return root;
}

void divideFacesByLocation(OctreeNode *node, std::list<DefaultMesh::FaceHandle> &facesList, const DefaultMesh &mesh)
{
    auto f_it = facesList.begin();
    //parent's inner triangles
    while(f_it != facesList.end())
    {
        switch (triangleRelativePosition(mesh, *f_it, node->min, node->size))
        {
            case INSIDE:
                node->innerFaces.push_back(*f_it);
                /*If the triangle is located inside the cell, we remove it from the cell's parent list*/
                f_it = facesList.erase(f_it);
                break;
            case CROSSING:
                node->crossingFaces.push_back(*f_it);
                /*If the triangle might cross the cell, we can't remove it from the cell's parent list
                 * because it might cross other cells as well*/
                ++f_it;
                break;
            default:
                ++f_it;
        }
    }
}

OctreeNode* ConstructOctreeNodesFromOpenMesh(OctreeNode *node, const DefaultMesh &mesh) {
    if (!node)
    {
        std::cout << "Trying to construct empty node" << std::endl;
        return nullptr;
    }

    if (node->parent != nullptr)
    {
        divideFacesByLocation(node, node->parent->innerFaces, mesh); //parent's inner triangles
        divideFacesByLocation(node, node->parent->crossingFaces, mesh); //parent's crossing triangles
    }
    else
    {   //initializes the parent list with all triangles
        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
        {
            node->innerFaces.push_back(*f_it);
        }
    }

    if (node->innerFaces.size() == 0)
    {   //Empty space, no triangles crossing or inside this cell
        if (node->crossingFaces.size() == 0)
        {
            //std::cout << "Empty space HERE!!" << std::endl;
            delete node;
            return nullptr;
        }
        else
        {
            if (node->parent->innerFaces.size() == 0)
                return ConstructLeafFromOpenMesh(node, mesh);
        }
    }

    if (node->height == 0)
    {
        //std::cout << "HEIGHT 0 ACTIVATED!!" << std::endl;
        return ConstructLeafFromOpenMesh(node, mesh);
    }

    const float childSize = node->size / 2;
    const int childHeight = node->height - 1;
    bool hasChildren = false;

    for (int i = 0; i < 8; i++)
    {
        OctreeNode* child = new OctreeNode;
        child->size = childSize;
        child->height = childHeight;
        child->min = node->min + (CHILD_MIN_OFFSETS[i] * childSize);
        child->type = Node_Internal;
        child->parent = node;

        node->children[i] = ConstructOctreeNodesFromOpenMesh(child, mesh);
        hasChildren |= (node->children[i] != nullptr);
    }

    if (!hasChildren)
    {
        delete node;
        return nullptr;
    }

    return node;
}

OctreeNode *ConstructLeafFromOpenMesh(OctreeNode *leaf, const DefaultMesh &mesh) {
    if (!leaf)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }
    //std::cout << "Leaf height: " << leaf->height << std::endl;
    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    svd::QefSolver featureQef;
    bool hasIntersection = false;
    int corners = 0;
    //vertices classification
    int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                       MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
    for (int i = 0; i < 12; i++) //for each edge
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

        //computes hash for edge
        std::string edgehash = hashedge(p1, p2);
        HermiteData edgedata;
        if (OctreeNode::edgepool.count(edgehash) != 0){
            // if edge data already exists, retrieve the data
            edgedata = OctreeNode::edgepool[edgehash];
            if (edgedata.hasIntersection()){
//                averageNormal += edgedata.normal;
//                qef.add(edgedata.intersection.x, edgedata.intersection.y, edgedata.intersection.z,
//                        edgedata.normal.x, edgedata.normal.y, edgedata.normal.z);
//                hasIntersection = true;
                vecsigns[c1] = OctreeNode::vertexpool[hashvertex(p1)];
                vecsigns[c2] = OctreeNode::vertexpool[hashvertex(p2)];
                if (vecsigns[c1] == MATERIAL_UNKNOWN || vecsigns[c2] == MATERIAL_UNKNOWN)
                {
                    std::cout << "SIGNAL INCONSISTENCY FOR VERTICES IN INTERSECTION" << std::endl;
                }
//                vecsigns[c1] = computeSideOfPoint(p1, edgedata.intersection, edgedata.normal);
//                vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
            }
            //continue;
        }

        vec3 intersection;
        std::vector<vec3> intersection_points;
        std::vector<vec3> normals;
        std::vector<vec3> face_normals;
        for (std::list<DefaultMesh::FaceHandle>::iterator face = leaf->crossingFaces.begin(); face != leaf->crossingFaces.end(); ++face)
        {
            auto fv_it = mesh.cfv_iter(*face);
            DefaultMesh::VertexHandle a = *fv_it;
            DefaultMesh::VertexHandle b = *(++fv_it);
            DefaultMesh::VertexHandle c = *(++fv_it);

            vec3 face_vertices[3] = {openmesh_to_glm(mesh.point(a)), openmesh_to_glm(mesh.point(b)), openmesh_to_glm(mesh.point(c))};
            Vertex vertices[3] = { face_vertices[0], face_vertices[1], face_vertices[2]};

            if (moller_triangle_intersection(p1, p2, vertices, intersection)) {
                //keeps the intersection here
                if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                    continue;
                }
                intersection_points.push_back(intersection);

                float u, v, w;
                barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);

                vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                normal_at_intersection =  glm::normalize(normal_at_intersection);
                normals.push_back(normal_at_intersection);
                vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                face_normals.push_back(face_normal);

                if (vecsigns[c1] != MATERIAL_SOLID/*vecsigns[c1] == MATERIAL_UNKNOWN*/){
                    vecsigns[c1] = computeSideOfPoint(p1, intersection, face_normal);
                }
                if (vecsigns[c2] != MATERIAL_SOLID/*vecsigns[c2] == MATERIAL_UNKNOWN*/){
                    vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
                }

                if (mesh.is_boundary(*face))
                {
                    for (auto fh_iter = mesh.cfh_iter(*face); fh_iter.is_valid(); ++fh_iter)
                    {
                        DefaultMesh::Normal faceNormal = mesh.normal(*face);
                        DefaultMesh::Point middle_point = (mesh.point(mesh.to_vertex_handle(*fh_iter)) + mesh.point(mesh.from_vertex_handle((*fh_iter))))/2;;
                        DefaultMesh::Point edge = mesh.point(mesh.to_vertex_handle(*fh_iter)) - mesh.point(mesh.from_vertex_handle((*fh_iter)));
                        DefaultMesh::Normal planeNormal;
                        if (mesh.is_boundary(*fh_iter))
                        {
                            planeNormal =  edge % faceNormal;
                            //planeNormal.normalize();
                            featureQef.add(middle_point[0], middle_point[1], middle_point[2], planeNormal[0], planeNormal[1], planeNormal[2]);
                        }
                        DefaultMesh::HalfedgeHandle opposite_halfedge = mesh.opposite_halfedge_handle(*fh_iter);
                        if (mesh.is_boundary(opposite_halfedge))
                        {
                            planeNormal = faceNormal % edge;
                            //planeNormal.normalize();
                            featureQef.add(middle_point[0], middle_point[1], middle_point[2], planeNormal[0], planeNormal[1], planeNormal[2]);
                        }
                    }
                }
            }
        }
        if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
            if (intersection_points.size()%2 == 0){
//                int nearindex = glm::distance(p1, intersection_points[0]) < glm::distance(p1, intersection_points[1]) ? 0 : 1;
//                vecsigns[c1] = computeSideOfPoint(p1, intersection_points[nearindex], face_normals[nearindex]);
//                // they are on the same side of the surface. We must ignore the intersection
//                vecsigns[c2] = vecsigns[c1];
                const float childSize = leaf->size / 2;
                const int childHeight = leaf->height - 1;
                std::cout << intersection_points.size() << " Child Height: " << childHeight << " Child Size: " << childSize << std:: endl;
                bool hasChildren = false;

                for (int i = 0; i < 8; i++)
                {
                    OctreeNode* child = new OctreeNode;
                    child->size = childSize;
                    child->height = childHeight;
                    child->min = leaf->min + (CHILD_MIN_OFFSETS[i] * childSize);
                    child->type = Node_Internal;
                    child->parent = leaf;

                    leaf->children[i] = ConstructOctreeNodesFromOpenMesh(child, mesh);
                    hasChildren |= (leaf->children[i] != nullptr);
                }

                if (!hasChildren)
                {
                    delete leaf;
                    return nullptr;
                }

                return leaf;
            }
            else{
                // they are on the opposite side of the surface
                int nearindex = glm::distance(p1, intersection_points[0]) < glm::distance(p1, intersection_points[1]) ? 0 : 1;
                nearindex = glm::distance(p1, intersection_points[nearindex]) < glm::distance(p1, intersection_points[2]) ? nearindex : 2;
                vecsigns[c1] = computeSideOfPoint(p1, intersection_points[nearindex], face_normals[nearindex]);
                // they are on the same side of the surface. We must ignore the intersection
                vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;

                std::cout << "CHEGA MAIS" << std::endl;
            }


        }
        // if we consider that an intersection happened.
        if ((intersection_points.size() > 0) && (vecsigns[c1] != vecsigns[c2]))
        {   // we'll consider only the first intersection for now
            vec3 &n = normals[0];
            vec3 &v = intersection_points[0];
            qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
            averageNormal += n;
            hasIntersection = true;
            edgedata.intersection = intersection_points[0];
            edgedata.normal = normals[0];
        }

        OctreeNode::edgepool[edgehash] = edgedata;
    }

    if (!hasIntersection)
    {   // voxel is full inside or outside the volume
        delete leaf;
        return nullptr;
    }
    updateSignsArray(vecsigns, 8);

    for (size_t i = 0; i < 8; i++)
    {   //encode the signs to the corners variable to save memory
        corners |= (vecsigns[i] << i);
        updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
    }


    std::ofstream interiorfile, exteriorfile;
    //gridfile.open("/home/hallpaz/Workspace/dual_contouring_experiments/grid_color_updated.ply", std::ios::app);
    interiorfile.open("../subproducts/interior_color_updated.ply", std::ios::app);
    exteriorfile.open("../subproducts/exterior_color_updated.ply", std::ios::app);
    for (size_t i = 0; i < 8; i++) {
        corners |= (vecsigns[i] << i);

        const vec3 cornerPos = leaf->min + CHILD_MIN_OFFSETS[i]*leaf->size;
        if (vecsigns[i] == MATERIAL_SOLID) {
            interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 128 << " " << 255 << std::endl;
        }

        if (vecsigns[i] == MATERIAL_AIR) {
            exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 64 << " " << 255 << " " << 64 << std::endl;
        }

        if (vecsigns[i] == MATERIAL_UNKNOWN) {
            interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
            exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
        }
    }
    interiorfile.close();
    exteriorfile.close();


    svd::Vec3 qefPosition;
    qef.setData(qef.getData()*0.4f + featureQef.getData()*0.6f);
    //qef.add(featureQef.getData());
    qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);

    OctreeDrawInfo* drawInfo = new OctreeDrawInfo;
    drawInfo->position = vec3(qefPosition.x, qefPosition.y, qefPosition.z);
    drawInfo->qef = qef.getData();

    const vec3 min = vec3(leaf->min);
    const vec3 max = vec3(leaf->min + vec3(leaf->size));
    if (drawInfo->position.x < min.x || drawInfo->position.x > max.x ||
        drawInfo->position.y < min.y || drawInfo->position.y > max.y ||
        drawInfo->position.z < min.z || drawInfo->position.z > max.z)
    {
        const auto& mp = qef.getMassPoint();
        drawInfo->position = vec3(mp.x, mp.y, mp.z);
    }

    drawInfo->averageNormal = glm::normalize(averageNormal);
    drawInfo->corners = corners;

    leaf->type = Node_Leaf;
    leaf->drawInfo = drawInfo;

    return leaf;
}
// -------------------------------------------------------------------------------

OctreeNode* SimplifyOctree(OctreeNode* node, const float threshold)
{
    if (!node)
    {
        //std::cout << "Empty node" << std::endl;
        return nullptr;
    }

    if (node->type != Node_Internal)
    {
        // can't simplify!
        return node;
    }
    //std::cout << "Simplifying at level: " << node->height << std::endl;
    svd::QefSolver qef;

    int signs[8] = { MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                    MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
    int midsign = MATERIAL_UNKNOWN;
    int edgeCount = 0;
    bool isCollapsible = true;

    for (int i = 0; i < 8; i++)
    {
        node->children[i] = SimplifyOctree(node->children[i], threshold);
        if (node->children[i])
        {
            OctreeNode* child = node->children[i];
            if (child->type == Node_Internal)
            {
                isCollapsible = false;
            }
            else
            {
                qef.add(child->drawInfo->qef);

                midsign = (child->drawInfo->corners >> (7 - i)) & 1;
                signs[i] = (child->drawInfo->corners >> i) & 1;

                edgeCount++;
            }
        }
    }

    if (!isCollapsible)
    {
        // at least one child is an internal node, can't collapse
        return node;
    }

    svd::Vec3 qefPosition;
    qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
    float error = qef.getError();

    // convert to glm vec3 for ease of use
    vec3 position(qefPosition.x, qefPosition.y, qefPosition.z);

    // at this point the masspoint will actually be a sum, so divide to make it the average
    if (error > threshold)
    {
        // this collapse breaches the threshold
        return node;
    }

    if ((position.x < node->min.x) || (position.x > (node->min.x + node->size)) ||
        (position.y < node->min.y) || (position.y > (node->min.y + node->size)) ||
        (position.z < node->min.z) || (position.z > (node->min.z + node->size)))
    {
        const auto& mp = qef.getMassPoint();
        position = vec3(mp.x, mp.y, mp.z);
    }

    // change the node from an internal node to a 'pseudo leaf' node
    OctreeDrawInfo* drawInfo = new OctreeDrawInfo;

    for (int i = 0; i < 8; i++)
    {
        if (signs[i] == MATERIAL_UNKNOWN)
        {
            // Undetermined, use centre sign instead
            drawInfo->corners |= (midsign << i);
        }
        else
        {
            drawInfo->corners |= (signs[i] << i);
        }
    }

    drawInfo->averageNormal = vec3(0.f);
    for (int i = 0; i < 8; i++)
    {
        if (node->children[i])
        {
            OctreeNode* child = node->children[i];
            if (child->type != Node_Internal /*child->type == Node_Pseudo || child->type == Node_Leaf*/)
            {
                drawInfo->averageNormal += child->drawInfo->averageNormal;
            }
        }
    }

    drawInfo->averageNormal = glm::normalize(drawInfo->averageNormal);
    drawInfo->position = position;
    drawInfo->qef = qef.getData();

    for (int i = 0; i < 8; i++)
    {
        DestroyOctree(node->children[i]);
        node->children[i] = nullptr;
    }

    node->type = Node_Pseudo;
    node->drawInfo = drawInfo;

    return node;
}

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

#include "../Constants.h"
#include "../Contouring.h"
#include "../density.h"
#include "octree.h"
#include "../Utils.h"

#include "../glm/ext.hpp"
/*

using glm::vec3;
using glm::ivec3;


// -------------------------------------------------------------------------------

std::unordered_map<std::string, int> OctreeNode::vertexpool;
std::unordered_map<std::string, HermiteData> OctreeNode::edgepool;

// ----------------------------------------------------------------------------

OctreeNode* ConstructLeaf(OctreeNode* leaf)
{

    if (!leaf || leaf->depth > 0)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }

    int corners = 0;
    for (int i = 0; i < 8; ++i)
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

    leaf->type = NODE_LEAF;
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

    if (node->depth == 0)
    {
        return construct_or_update_leaf(node);
    }

    const float childSize = node->size / 2;
    const int childHeight = node->depth - 1;
    bool hasChildren = false;

    for (int i = 0; i < 8; i++)
    {
        OctreeNode* child = new OctreeNode;
        child->size = childSize;
        child->depth = childHeight;
        child->min = node->min + (CHILD_MIN_OFFSETS[i] * childSize);
        child->type = NODE_INTERNAL;

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
    root->depth = height;
    root->type = NODE_INTERNAL;

    std::cout << "Octree.cpp: will construct nodes" << std::endl;
    ConstructOctreeNodes(root);
    std::cout << "Octree.cpp: will simplify nodes" << std::endl;
    root = SimplifyOctree(root, threshold);
    std::cout << "Octree.cpp: did simplify nodes" << std::endl;

    return root;
}

// -------------------------------------------------------------------------------

OctreeNode::~OctreeNode()
{
    delete drawInfo;
    for (int i = 0; i < 8; i++)
    {
        delete children[i];
    }
}

// -------------------------------------------------------------------------------

OctreeNode* BuildOctreeFromOpenMesh(const glm::vec3 &min, const float size, const int height, const DefaultMesh &mesh)
{
    trace("Octree.cpp: BuildOctreeFromOpenMesh");
    OctreeNode* root = new OctreeNode(NODE_INTERNAL, min, size, height);
    ConstructOctreeNodesFromOpenMesh(root, mesh);
    return root;
}

OctreeNode* ConstructOctreeNodesFromOpenMesh(OctreeNode *node, const DefaultMesh &mesh)
{
    if (!node)
    {
        std::cout << "Trying to construct empty node" << std::endl;
        return nullptr;
    }
    select_inner_crossing_faces(node, mesh);
    if (node->innerFaces.empty() && node->crossingFaces.empty())
    {   //Empty space, no triangles crossing or inside this cell
        delete node;
        return nullptr;
    }
    if ((node->parent && node->parent->innerFaces.empty()) || node->depth == 0)
    {
        return ConstructLeafFromOpenMesh(node, mesh);
    }

    bool has_children = construct_or_update_children(node, mesh);
    if (has_children)
    {
        //clear memory used for inner and crossing faces
        return clean_node(node);
    }

    delete node;
    return nullptr;
}

OctreeNode *ConstructLeafFromOpenMesh(OctreeNode *leaf, const DefaultMesh &mesh) {
    if (!leaf)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }

    for (int j = 0; j < 8; ++j) {
        OctreeNode::vertexpool[hashvertex(vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[j]))] = MATERIAL_UNKNOWN;
    }

    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    bool hasIntersection = false;
    int corners = 0;
    //vertices classification
    int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                       MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
    for (int i = 0; i < 12; ++i) //for each edge
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
                vecsigns[c1] = OctreeNode::vertexpool[hashvertex(p1)];
                vecsigns[c2] = OctreeNode::vertexpool[hashvertex(p2)];
                if (vecsigns[c1] == MATERIAL_UNKNOWN || vecsigns[c2] == MATERIAL_UNKNOWN)
                {
                    std::cout << "SIGNAL INCONSISTENCY FOR VERTICES IN INTERSECTION" << std::endl;
                }
            }
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
                hasIntersection = true;
            }
        }
    }

    if (!hasIntersection)
    {   // voxel is full inside or outside the volume
        delete leaf;
        return nullptr;
    }
    for (size_t i = 0; i < 8; i++)
    {   //encode the signs to the corners variable to save memory
        corners |= (vecsigns[i] << i);
        updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
    }

    // DEBUG --------------------------------------------------------------------
    std::ofstream interiorfile, exteriorfile;
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
    // DEBUG --------------------------------------------------------------------//

    OctreeDrawInfo* drawInfo = new OctreeDrawInfo;
    drawInfo->qef = qef.getData();
    drawInfo->averageNormal += averageNormal;
    drawInfo->corners = corners;
    leaf->drawInfo = drawInfo;

    leaf->type = NODE_LEAF;

    //return clean_node(leaf);
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

    if (node->type != NODE_INTERNAL)
    {
        // can't simplify!
        return node;
    }
    //std::cout << "Simplifying at level: " << node->depth << std::endl;
    svd::QefSolver qef;

    int signs[8] = { MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                    MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN };
    int midsign = MATERIAL_UNKNOWN;
    int edgeCount = 0;
    bool isCollapsible = true;

    for (int i = 0; i < 8; ++i)
    {
        node->children[i] = SimplifyOctree(node->children[i], threshold);
        if (node->children[i])
        {
            OctreeNode* child = node->children[i];
            if (child->type == NODE_INTERNAL)
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
    for (int i = 0; i < 8; ++i)
    {
        if (node->children[i])
        {
            OctreeNode* child = node->children[i];
            if (child->type != NODE_INTERNAL)
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
        delete node->children[i];
        node->children[i] = nullptr;
    }

    node->type = NODE_PSEUDO;
    node->drawInfo = drawInfo;

    return node;
}
*/
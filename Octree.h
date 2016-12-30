//
// Created by Hallison da Paz on 18/11/2016.
//
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
#ifndef DUAL_CONTOURING_SANDBOX_OCTREE_H
#define DUAL_CONTOURING_SANDBOX_OCTREE_H

#include <vector>
#include <unordered_map>
#include "glm/vec3.hpp"

#include "Constants.h"
#include "DataStructures.h"


class OctreeNode
{
    friend class Octree;
public:

    OctreeNode(OctreeNodeType type, glm::vec3 min, float size, int height, OctreeNode* parent = nullptr)
            : type(type), min(min), size(size), depth(height), drawInfo(nullptr), parent(parent)
    {
        for (int i = 0; i < NUM_CHILDREN; ++i)
        {
            children[i] = nullptr;
        }
    }

    OctreeNode() : type(NODE_NONE), min(0, 0, 0), size(0), depth(0), drawInfo(nullptr), parent(nullptr)
    {
        for (int i = 0; i < 8; ++i)
        {
            children[i] = nullptr;
        }
    }

    ~OctreeNode()
    {
        for (int i = 0; i < NUM_CHILDREN; ++i)
        {
            delete children[i];
        }
        delete drawInfo;
        //delete meshInfo;
    }

    bool construct_or_update_children(unsigned int max_depth, const DefaultMesh &mesh);

    void clean()
    {
        crossingFaces.clear();
        innerFaces.clear();
        /*delete meshInfo;
        meshInfo = nullptr;*/
    }

    bool inline crossingEmpty()
    {
        /*if (meshInfo == nullptr)
            return true;
        return meshInfo->crossingFaces.empty();*/
        return crossingFaces.empty();
    }

    bool inline innerEmpty()
    {
        /*if (meshInfo == nullptr)
            return true;
        return meshInfo->innerFaces.empty();*/
        return innerFaces.empty();
    }

    glm::vec3 get_vertex(int index)
    {
        assert(index >= 0 && index <= NUM_CHILDREN);
        return min + size * CHILD_MIN_OFFSETS[index];
    }

    glm::vec3 get_max()
    {
        return min + size*glm::vec3(1.0f, 1.0f, 1.0f);
    }

    OctreeNodeType	type;
    glm::vec3		min;
    Real			size;
    unsigned int    depth;
    OctreeNode*		children[8];
    OctreeNode*     parent;
    OctreeDrawInfo*	drawInfo;

    bool irregular = false; //irregular intersection
    int num_intersections = 0;
    bool dangerous_intersection = false;
    //OctreeMeshInfo* meshInfo;
    std::list<DefaultMesh::FaceHandle> innerFaces;
    std::list<DefaultMesh::FaceHandle> crossingFaces;

};


class Octree {
public:

    OctreeNode *root;
    Octree(glm::vec3 min, Real size, unsigned int max_depth, DefaultMesh &mesh, glm::vec3 cam_origin);
    void classify_leaves_vertices(glm::vec3 cam_origin, OctreeNode* node, DefaultMesh &mesh);
    static OctreeNode *BuildMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh);
    static OctreeNode *UpdateMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh);
    static OctreeNode *construct_or_update_leaf(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh);
    static void classify_leaves_vertices(OctreeNode* node);

    static std::unordered_map<std::string, int> leafvertexpool;
    static OctreeNode* SimplifyOctree(OctreeNode* node, const float threshold);


#ifdef DEBUG
    static int unoptimized_points;
    static int divergence;
    static int ambiguous_vertices;
    static int irregular_cells;
#endif

};


#endif //DUAL_CONTOURING_SANDBOX_OCTREE_H
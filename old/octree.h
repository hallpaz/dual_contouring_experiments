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


#ifndef		HAS_OCTREE_H_BEEN_INCLUDED
#define		HAS_OCTREE_H_BEEN_INCLUDED

#include <vector>
#include <list>
#include <unordered_map>

#include "qef.h"

#include "../glm/glm.hpp"
#include "../DataStructures.h"
#include "../Constants.h"

// ----------------------------------------------------------------------------

enum OctreeNodeType
{
    NODE_NONE,
    NODE_INTERNAL,
    NODE_PSEUDO,
    NODE_LEAF,
};

// ----------------------------------------------------------------------------

struct OctreeDrawInfo
{
    OctreeDrawInfo()
            : index(-1)
            , corners(255)
    {
    }

    int	index;
    int	corners;
    glm::vec3 position;
    glm::vec3 averageNormal;
    svd::QefData qef;
};

/*struct HermiteData
{
    glm::vec3 intersection;
    glm::vec3 normal;
    bool hasIntersection()
    {
        if (normal.x || normal.y || normal.z)
        {
            return true;
        }
        return false;
    }
};

// ----------------------------------------------------------------------------

class OctreeNode
{
public:

    OctreeNode(OctreeNodeType type, glm::vec3 min, float size, int height, OctreeNode* parent = nullptr)
            : type(type)
            , min(min)
            , size(size)
            , depth(height)
            , drawInfo(nullptr)
            , parent(parent)
    {
        for (int i = 0; i < 8; ++i)
        {
            children[i] = nullptr;
        }
    }

    OctreeNode()
            : type(NODE_NONE)
            , min(0, 0, 0)
            , size(0)
            , depth(0)
            , drawInfo(nullptr)
            , parent(nullptr)
    {
        for (int i = 0; i < 8; ++i)
        {
            children[i] = nullptr;
        }
    }

    OctreeNode(const OctreeNodeType _type)
            : type(_type)
            , min(0, 0, 0)
            , size(0)
            , depth(0)
            , drawInfo(nullptr)
            , parent(nullptr)
    {
        for (int i = 0; i < 8; i++)
        {
            children[i] = nullptr;
        }
    }

    ~OctreeNode();

    OctreeNodeType	type;
    glm::vec3		min;
    float			size;
    int             depth;
    OctreeNode*		children[8];
    OctreeNode*     parent;
    OctreeDrawInfo*	drawInfo;

    static std::unordered_map<std::string, int> vertexpool;
    static std::unordered_map<std::string, HermiteData> edgepool;

// ---------------------------------------------------------------------------- OPENMESH

    std::list<DefaultMesh::FaceHandle> innerFaces;
    std::list<DefaultMesh::FaceHandle> crossingFaces;

};

// ----------------------------------------------------------------------------

OctreeNode* BuildOctree(const glm::vec3& min, const float size, const int height, const float threshold);
OctreeNode* BuildOctreeFromOpenMesh(const glm::vec3& min, const float size, const int height, const DefaultMesh &mesh);

/* OPENMESH
OctreeNode *ConstructOctreeNodesFromOpenMesh(OctreeNode *pNode, const DefaultMesh &mesh);
OctreeNode *ConstructLeafFromOpenMesh(OctreeNode *node, const DefaultMesh &mesh);
OctreeNode* SimplifyOctree(OctreeNode* node, const float threshold);

// ----------------------------------------------------------------------------

inline bool construct_children(OctreeNode* node, const DefaultMesh &mesh)
{
    const float childSize = node->size / 2;
    const int childHeight = node->depth - 1;
    bool hasChildren = false;
    for (int i = 0; i < 8; i++)
    {
        OctreeNode* child = new OctreeNode(NODE_INTERNAL,
                                           node->min + (CHILD_MIN_OFFSETS[i] * childSize),
                                           childSize,
                                           childHeight,
                                           node);
        node->children[i] = ConstructOctreeNodesFromOpenMesh(child, mesh);
        hasChildren |= (node->children[i] != nullptr);
    }
    return hasChildren;
}

inline OctreeNode* clean_node(OctreeNode* node)
{
    node->innerFaces.clear();
    node->crossingFaces.clear();
    return node;
}*/

// ----------------------------------------------------------------------------

#endif	// HAS_OCTREE_H_BEEN_INCLUDED
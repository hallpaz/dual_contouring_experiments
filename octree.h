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

#include "qef.h"

#include "glm/glm.hpp"
#include "dataStructures.h"
using glm::vec3;
using glm::ivec3;

// ----------------------------------------------------------------------------

enum OctreeNodeType
{
    Node_None,
    Node_Internal,
    Node_Pseudo,
    Node_Leaf,
};

enum RelativePosition
{
    INSIDE,
    CROSSING,
    OUTSIDE,
};

// ----------------------------------------------------------------------------

struct OctreeDrawInfo
{
    OctreeDrawInfo()
            : index(-1)
            , corners(0)
    {
    }

    int				index;
    int				corners;
    vec3			position;
    vec3			averageNormal;
    svd::QefData	qef;
};

// ----------------------------------------------------------------------------

class OctreeNode
{
public:

    OctreeNode()
            : type(Node_None)
            , min(0, 0, 0)
            , size(0)
            , height(0)
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
            , height(0)
            , drawInfo(nullptr)
            , parent(nullptr)
    {
        for (int i = 0; i < 8; i++)
        {
            children[i] = nullptr;
        }
    }

    OctreeNodeType	type;
    vec3			min;
    float			size;
    int             height;
    OctreeNode*		children[8];
    OctreeNode*     parent;
    OctreeDrawInfo*	drawInfo;
    std::list<Triangle > innerTriangles;
    std::list<Triangle > crossingTriangles;

    RelativePosition vertexRelativePosition(const Vertex &vertex);
    RelativePosition triangleRelativePosition(const Vertex &a, const Vertex& b, const Vertex& c);
};

// ----------------------------------------------------------------------------

OctreeNode* BuildOctree(const vec3& min, const float size, const int height, const float threshold);
void DestroyOctree(OctreeNode* node);
void GenerateMeshFromOctree(OctreeNode* node, VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer);

OctreeNode* BuildOctreeFromMesh(const vec3& min, const float size, const int height, const float threshold,
                                VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer);



// ----------------------------------------------------------------------------

#endif	// HAS_OCTREE_H_BEEN_INCLUDED
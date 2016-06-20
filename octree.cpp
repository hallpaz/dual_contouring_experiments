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

#include	"octree.h"
#include	"density.h"
#include <fstream>
#include <sstream>
#include <iomanip>

// ----------------------------------------------------------------------------

const int MATERIAL_AIR = 0;
const int MATERIAL_SOLID = 1;
const int MATERIAL_UNKNOWN = -1;

const float QEF_ERROR = 1e-6f;
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

const vec3 CHILD_MIN_OFFSETS[] =
        {
                // needs to match the vertMap from Dual Contouring impl
                vec3( 0, 0, 0 ),
                vec3( 0, 0, 1 ),
                vec3( 0, 1, 0 ),
                vec3( 0, 1, 1 ),
                vec3( 1, 0, 0 ),
                vec3( 1, 0, 1 ),
                vec3( 1, 1, 0 ),
                vec3( 1, 1, 1 ),
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

OctreeNode* ConstructOctreeNodesFromMesh(OctreeNode *pNode, const VertexBuffer &vector, const IndexBuffer &buffer);

OctreeNode *ConstructLeafFromMesh(OctreeNode *node, const VertexBuffer &vertexBuffer, const IndexBuffer &indexBuffer);

std::unordered_map<std::string, int> OctreeNode::vertexpool;
std::unordered_map<std::string, HermiteData> OctreeNode::edgepool;

std::string hashvertex(const vec3& vertex){
    int precision = 3;
    std::stringstream rep;
    rep << vertex.x << std::setprecision (precision) << std::fixed
        << vertex.y << std::setprecision (precision) << std::fixed
        << vertex.z << std::setprecision (precision) << std::fixed;
    return rep.str();
}

std::string hashedge(const vec3& v1, const vec3& v2){
    if (v1.x < v2.x){
        return (hashvertex(v1) + hashvertex(v2));
    }
    return (hashvertex(v2) + hashvertex(v1));
}

OctreeNode* SimplifyOctree(OctreeNode* node, float threshold)
{
    if (!node)
    {
        //std::cout << "Empty node" << std::endl;
        return NULL;
    }

    if (node->type != Node_Internal)
    {
        // can't simplify!
        return node;
    }
    //std::cout << "Simplifying at level: " << node->height << std::endl;
    svd::QefSolver qef;

    int signs[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
    int midsign = -1;
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
        if (signs[i] == -1)
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
            if (child->type == Node_Pseudo ||
                child->type == Node_Leaf)
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

// ----------------------------------------------------------------------------

void GenerateVertexIndices(OctreeNode* node, VertexBuffer& vertexBuffer)
{
    if (!node)
    {
        return;
    }

    if (node->type != Node_Leaf)
    {
        for (int i = 0; i < 8; i++)
        {
            GenerateVertexIndices(node->children[i], vertexBuffer);
        }
    }

    if (node->type != Node_Internal)
    {
        OctreeDrawInfo* d = node->drawInfo;
        if (!d)
        {
            std::cout << "Error! Could not add vertex!" << std::endl;
            exit(EXIT_FAILURE);
        }

        d->index = vertexBuffer.size();
        //vertexBuffer.push_back(d->position);
        vertexBuffer.push_back(Vertex(d->position));
    }
}

// ----------------------------------------------------------------------------

void ContourProcessEdge(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer)
{
    float minSize = 1000000.f;		// arbitrary big number
    int minIndex = 0;
    int indices[4] = { -1, -1, -1, -1 };
    bool flip = false;
    bool signChange[4] = { false, false, false, false };

    for (int i = 0; i < 4; i++)
    {
        const int edge = processEdgeMask[dir][i];
        const int c1 = edgevmap[edge][0];
        const int c2 = edgevmap[edge][1];

        const int m1 = (node[i]->drawInfo->corners >> c1) & 1;
        const int m2 = (node[i]->drawInfo->corners >> c2) & 1;

        if (node[i]->size < minSize)
        {
            minSize = node[i]->size;
            minIndex = i;
            flip = m1 != MATERIAL_AIR;
        }

        indices[i] = node[i]->drawInfo->index;

        signChange[i] =
                (m1 == MATERIAL_AIR && m2 != MATERIAL_AIR) ||
                (m1 != MATERIAL_AIR && m2 == MATERIAL_AIR);
    }

    if (signChange[minIndex])
    {
        if (!flip)
        {
            /*indexBuffer.push_back(indices[0]);
            indexBuffer.push_back(indices[1]);
            indexBuffer.push_back(indices[3]);*/
            indexBuffer.push_back(Triangle{indices[0], indices[1], indices[3]});

            /*indexBuffer.push_back(indices[0]);
            indexBuffer.push_back(indices[3]);
            indexBuffer.push_back(indices[2]);*/
            indexBuffer.push_back(Triangle{indices[0], indices[3], indices[2]});
        }
        else
        {
            /*indexBuffer.push_back(indices[0]);
            indexBuffer.push_back(indices[3]);
            indexBuffer.push_back(indices[1]);*/
            indexBuffer.push_back(Triangle{indices[0], indices[3], indices[1]});

            /*indexBuffer.push_back(indices[0]);
            indexBuffer.push_back(indices[2]);
            indexBuffer.push_back(indices[3]);*/
            indexBuffer.push_back(Triangle{indices[0], indices[2], indices[3]});
        }
    }
}

// ----------------------------------------------------------------------------

void ContourEdgeProc(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer)
{
    if (!node[0] || !node[1] || !node[2] || !node[3])
    {
        return;
    }

    if (node[0]->type != Node_Internal &&
        node[1]->type != Node_Internal &&
        node[2]->type != Node_Internal &&
        node[3]->type != Node_Internal)
    {
        ContourProcessEdge(node, dir, indexBuffer);
    }
    else
    {
        for (int i = 0; i < 2; i++)
        {
            OctreeNode* edgeNodes[4];
            const int c[4] =
                    {
                            edgeProcEdgeMask[dir][i][0],
                            edgeProcEdgeMask[dir][i][1],
                            edgeProcEdgeMask[dir][i][2],
                            edgeProcEdgeMask[dir][i][3],
                    };

            for (int j = 0; j < 4; j++)
            {
                if (node[j]->type == Node_Leaf || node[j]->type == Node_Pseudo)
                {
                    edgeNodes[j] = node[j];
                }
                else
                {
                    edgeNodes[j] = node[j]->children[c[j]];
                }
            }

            ContourEdgeProc(edgeNodes, edgeProcEdgeMask[dir][i][4], indexBuffer);
        }
    }
}

// ----------------------------------------------------------------------------

void ContourFaceProc(OctreeNode* node[2], int dir, IndexBuffer& indexBuffer)
{
    if (!node[0] || !node[1])
    {
        return;
    }

    if (node[0]->type == Node_Internal ||
        node[1]->type == Node_Internal)
    {
        for (int i = 0; i < 4; i++)
        {
            OctreeNode* faceNodes[2];
            const int c[2] =
                    {
                            faceProcFaceMask[dir][i][0],
                            faceProcFaceMask[dir][i][1],
                    };

            for (int j = 0; j < 2; j++)
            {
                if (node[j]->type != Node_Internal)
                {
                    faceNodes[j] = node[j];
                }
                else
                {
                    faceNodes[j] = node[j]->children[c[j]];
                }
            }

            ContourFaceProc(faceNodes, faceProcFaceMask[dir][i][2], indexBuffer);
        }

        const int orders[2][4] =
                {
                        { 0, 0, 1, 1 },
                        { 0, 1, 0, 1 },
                };
        for (int i = 0; i < 4; i++)
        {
            OctreeNode* edgeNodes[4];
            const int c[4] =
                    {
                            faceProcEdgeMask[dir][i][1],
                            faceProcEdgeMask[dir][i][2],
                            faceProcEdgeMask[dir][i][3],
                            faceProcEdgeMask[dir][i][4],
                    };

            const int* order = orders[faceProcEdgeMask[dir][i][0]];
            for (int j = 0; j < 4; j++)
            {
                if (node[order[j]]->type == Node_Leaf ||
                    node[order[j]]->type == Node_Pseudo)
                {
                    edgeNodes[j] = node[order[j]];
                }
                else
                {
                    edgeNodes[j] = node[order[j]]->children[c[j]];
                }
            }

            ContourEdgeProc(edgeNodes, faceProcEdgeMask[dir][i][5], indexBuffer);
        }
    }
}

// ----------------------------------------------------------------------------

void ContourCellProc(OctreeNode* node, IndexBuffer& indexBuffer)
{
    if (node == NULL)
    {
        return;
    }

    if (node->type == Node_Internal)
    {
        for (int i = 0; i < 8; i++)
        {
            ContourCellProc(node->children[i], indexBuffer);
        }

        for (int i = 0; i < 12; i++)
        {
            OctreeNode* faceNodes[2];
            const int c[2] = { cellProcFaceMask[i][0], cellProcFaceMask[i][1] };

            faceNodes[0] = node->children[c[0]];
            faceNodes[1] = node->children[c[1]];

            ContourFaceProc(faceNodes, cellProcFaceMask[i][2], indexBuffer);
        }

        for (int i = 0; i < 6; i++)
        {
            OctreeNode* edgeNodes[4];
            const int c[4] =
                    {
                            cellProcEdgeMask[i][0],
                            cellProcEdgeMask[i][1],
                            cellProcEdgeMask[i][2],
                            cellProcEdgeMask[i][3],
                    };

            for (int j = 0; j < 4; j++)
            {
                edgeNodes[j] = node->children[c[j]];
            }

            ContourEdgeProc(edgeNodes, cellProcEdgeMask[i][4], indexBuffer);
        }
    }
}

// ----------------------------------------------------------------------------
//TODO: give a better solution to the steps MAGIC NUMBER
/*vec3 ApproximateZeroCrossingPosition(const vec3& p0, const vec3& p1, float size)
{
    // approximate the zero crossing by finding the min value along the edge
    float minValue = 100000.f;
    float t = 0.f;
    float currentT = 0.f;
    const int steps = 8;
    const float increment = size / (float)steps;
    while (currentT <= size)
    {
        const vec3 p = p0 + ((p1 - p0) * currentT);
        const float density = glm::abs(Density_Func(p));
        if (density < minValue)
        {
            minValue = density;
            t = currentT;
        }

        currentT += increment;
    }

    return p0 + ((p1 - p0) * t);
}*/

vec3 ApproximateZeroCrossingPosition(const vec3& p0, const vec3& p1, float size)
{
    const float tolerance = 0.00001;
    const int maxIterations = 21;
    float tbegin = 0.0;
    float tend = 1.0;
    float t = (tbegin + tend)/2;
    int currentIt = 0;
    float middle = Density_Func(p0 + ((p1 - p0) * t));
    float begin = Density_Func(p0 + ((p1 - p0) * tbegin));
    float end = Density_Func(p0 + ((p1 - p0) * tend));
    //std::cout << "Begin: " << begin << " Middle: " << middle << " End: " << end << std::endl;
    while (currentIt < maxIterations)
    {

        if (fabs(middle) < fabs(tolerance)){
            break;
        }
        //double check
        if (begin * end > 0){
            std::cout << "Weird, but this is not a bipolar edge! it: " << currentIt << std::endl;
        }

        //check the side of the middle point
        if(middle*begin > 0){ //same side of begin
            tbegin = t;
        }
        else { //same side of end
            tend = t;
        }
        t = (tbegin + tend)/2;
        middle = Density_Func(p0 + ((p1 - p0) * t));
        begin = Density_Func(p0 + ((p1 - p0) * tbegin));
        end = Density_Func(p0 + ((p1 - p0) * tend));
        //std::cout << "Begin: " << begin << " Middle: " << middle << " End: " << end << std::endl;
    }
    std::cout << "Approximation: " << middle << std::endl;
    //exit(1);
    return (p0 + ((p1 - p0) * t));
}

// ----------------------------------------------------------------------------

vec3 CalculateSurfaceNormal(const vec3& p)
{
    const float H = 0.001f;
    const float dx = Density_Func(p + vec3(H, 0.f, 0.f)) - Density_Func(p - vec3(H, 0.f, 0.f));
    const float dy = Density_Func(p + vec3(0.f, H, 0.f)) - Density_Func(p - vec3(0.f, H, 0.f));
    const float dz = Density_Func(p + vec3(0.f, 0.f, H)) - Density_Func(p - vec3(0.f, 0.f, H));

    //return glm::normalize(vec3(dx, dy, dz));
    return glm::normalize(vec3(2*p.x, 2*p.y, 2*p.z));
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
        const vec3 p = ApproximateZeroCrossingPosition(p1, p2, leaf->size);
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

    /*if (node->size == 1)
    {
        return ConstructLeaf(node);
    }*/

    if (node->height == 0)
    {
        //std::cout << "Will construct Leaf" << std::endl;
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
        //std::cout << hasChildren << std::endl;
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

// ----------------------------------------------------------------------------

void GenerateMeshFromOctree(OctreeNode* node, VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer)
{
    if (!node)
    {
        return;
    }

    vertexBuffer.clear();
    indexBuffer.clear();

    GenerateVertexIndices(node, vertexBuffer);
    ContourCellProc(node, indexBuffer);
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
OctreeNode* BuildOctreeFromMesh(const vec3& min, const float size, const int height, const float threshold,
                                VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer){


    OctreeNode* root = new OctreeNode;
    root->min = min;
    root->size = size;
    root->height = height;
    root->type = Node_Internal;

    std::cout << "Octree.cpp: will construct nodes" << std::endl;
    ConstructOctreeNodesFromMesh(root, vertexBuffer, indexBuffer);
    /*std::cout << "Octree.cpp: will simplify nodes" << std::endl;
    root = SimplifyOctree(root, threshold);
    std::cout << "Octree.cpp: did simplify nodes" << std::endl;*/

    return root;

}

OctreeNode* ConstructOctreeNodesFromMesh(OctreeNode *node, const VertexBuffer &vector, const IndexBuffer &buffer) {
    if (!node)
    {
        std::cout << "Trying to construct empty node" << std::endl;
        return nullptr;
    }

    if (node->parent != nullptr){
        OctreeNode *parent = node->parent;
        auto titerator = parent->innerTriangles.begin();
        //parent's inner triangles
        while(titerator != parent->innerTriangles.end())
        {
            switch (node->triangleRelativePosition(vector[(*titerator).a], vector[(*titerator).b], vector[(*titerator).c]))
            {
                case INSIDE:
                    node->innerTriangles.push_back(*titerator);
                    titerator = parent->innerTriangles.erase(titerator);
                    break;
                case CROSSING:
                    node->crossingTriangles.push_back(*titerator);
                    //titerator = parent->innerTriangles.erase(titerator);
                    ++titerator;
                    break;
                default:
                    ++titerator;
            }
        }
        //parent's crossing triangles
        titerator = parent->crossingTriangles.begin();
        while(titerator != parent->crossingTriangles.end())
        {
            switch (node->triangleRelativePosition(vector[(*titerator).a], vector[(*titerator).b], vector[(*titerator).c]))
            {
                case INSIDE:
                    node->innerTriangles.push_back(*titerator);
                    titerator = parent->crossingTriangles.erase(titerator);
                    break;
                case CROSSING:
                    node->crossingTriangles.push_back(*titerator);
                    //titerator = parent->crossingTriangles.erase(titerator);
                    break;
                default:
                    ;
            }
            ++titerator;
        }
    }else {
        //initializes the parent list with all triangles
        for (auto titerator = buffer.begin(); titerator != buffer.end() ; ++titerator) {
            RelativePosition tposition = node->triangleRelativePosition(vector[(*titerator).a], vector[(*titerator).b], vector[(*titerator).c]);
           //std::cout << tposition << std::endl;
            //exit(654);
            switch (tposition)
            {
                case INSIDE:
                    node->innerTriangles.push_back(*titerator);
                    break;
                case CROSSING:
                    node->crossingTriangles.push_back(*titerator);
                    break;
                default:
                    break;
            }
        }
    }

//    std::cout << "Height: " << node->height << " Inner triangles: " << node->innerTriangles.size() << " Crossing Triangles: " << node->crossingTriangles.size() << std::endl;
    if (node->height == 0)
    {
        //std::cout << "Will construct Leaf" << std::endl;
        return ConstructLeafFromMesh(node, vector, buffer);
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

        node->children[i] = ConstructOctreeNodesFromMesh(child, vector, buffer);
        hasChildren |= (node->children[i] != nullptr);
        //std::cout << hasChildren << std::endl;
    }

    if (!hasChildren)
    {
        delete node;
        return nullptr;
    }

    return node;
}

bool moller_triangle_intersection(vec3 v1, vec3 v2, Vertex* triangle_vertices, vec3& intersection_point)
{
    //std::cout << "moller" << std::endl;
    float EPSILON = 0.000001;

    vec3 e1 = triangle_vertices[1].position - triangle_vertices[0].position;
    vec3 e2 = triangle_vertices[2].position - triangle_vertices[0].position;

    vec3 D = v2 - v1;
    auto P = glm::cross(D, e2);

    float det = glm::dot(e1, P);
    if (glm::abs(det) < EPSILON){
        return false;
    }


    float inv_det = 1.0/det;
    vec3 T = v1 - triangle_vertices[0].position;
    float u = glm::dot(T, P) * inv_det;

    if(u < 0.0f or u > 1.0f){
        return false;
    }


    vec3 Q = glm::cross(T, e1);
    float v = glm::dot(D, Q) * inv_det;
    if(v < 0.0f or (u + v)  > 1.0f) {
        return false;
    }

    float t = glm::dot(e2, Q) * inv_det;

    if(t > 0.0 && t < 1.0){
        //std::cout << "Valor de t: " << t << std::endl;
        intersection_point = v1 + (v2-v1)*t;
        return true;
    }/*else {
        std::cout << "Valor de t: " << t << std::endl;
        std::cout << "Fora do alcance do segmento" << std::endl;
    }*/


    return false;
}


vec3 CalculateMeshNormal(Vertex* triangleVertices)
{
    vec3 v1 = triangleVertices[1].position - triangleVertices[0].position;
    vec3 v2 = triangleVertices[2].position - triangleVertices[1].position;
    auto normal = glm::normalize(glm::cross(v1, v2));

    return normal;
}

OctreeNode *ConstructLeafFromMesh(OctreeNode *leaf, const VertexBuffer &vertexBuffer, const IndexBuffer &indexBuffer) {
    if (!leaf || leaf->height != 0)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }

    /*if (leaf->innerTriangles.size() > 0){
        std::cout << "Leaf size: " << leaf->size << "Inner Triangles: " << leaf->innerTriangles.size() << std::endl;
    }*/

    // otherwise the voxel contains the surface, so find the edge intersections
    int edgecount = 0;
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    bool hasIntersection = false;
    int corners = 0;
    int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                       MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};


    for (int i = 0; i < 12; i++)
    {
        //for (std::vector<Triangle>::const_iterator face = indexBuffer.begin(); face < indexBuffer.end(); ++face) {
        int numIntersectionsOnEdge = 0;

        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

        //computes hash for edge
        std::string edgehash = hashedge(p1, p2);
        HermiteData edgedata;
        if (OctreeNode::edgepool.count(edgehash) != 0){
            // if exists, retrieve the data
            edgedata = OctreeNode::edgepool[edgehash];
            if (edgedata.hasIntersection){
                averageNormal += edgedata.normal;
                qef.add(edgedata.intersection.x, edgedata.intersection.y, edgedata.intersection.z,
                        edgedata.normal.x, edgedata.normal.y, edgedata.normal.z);
                hasIntersection = true;

                vecsigns[c1] = OctreeNode::vertexpool[hashvertex(p1)];
                vecsigns[c2] = OctreeNode::vertexpool[hashvertex(p2)];
            }
            continue;
        }

        vec3 intersection(0.f);
        std::vector<vec3> intersection_points;
        std::vector<vec3> normals;

        for (std::list<Triangle>::iterator face = leaf->crossingTriangles.begin(); face != leaf->crossingTriangles.end(); ++face) {

            Vertex vertices[3] = { vertexBuffer[(*face).a], vertexBuffer[(*face).b], vertexBuffer[(*face).c]  };

            if (moller_triangle_intersection(p1, p2, vertices, intersection)) {
                //keeps the intersection here
                if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < 0.000001f)){
                    continue;
                }
                intersection_points.push_back(intersection);

                ++numIntersectionsOnEdge;
                const vec3 n = CalculateMeshNormal(vertices);
                normals.push_back(n);
                //qef.add(intersection.x, intersection.y, intersection.z, n.x, n.y, n.z);

                //averageNormal += n;

                //edgecount++;
                //std::cout << "edgecount: " << edgecount << " e1: " << c1 << " e2: " << c2 << std::endl;
                //hasIntersection = true;

                /*edgedata.hasIntersection = true;
                edgedata.intersection = intersection;
                edgedata.normal = n;*/

                const int sign1 = glm::dot((p1 - intersection), n) < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
                const int sign2 = sign1 == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;

                vecsigns[c1] = sign1;
                vecsigns[c2] = sign2;

                std::ofstream outfile;
                outfile.open("../../intersection_color.ply", std::ios::app);
                if (i < 4) { // x axis
                    outfile << intersection.x << " " << intersection.y << " " << intersection.z << " " << 255 << " " <<
                    255 << " " << 0 << std::endl;
                }
                else {
                    if (i < 8) { // y axis
                        outfile << intersection.x << " " << intersection.y << " " << intersection.z << " " << 0 <<
                        " " << 255 << " " << 255 << std::endl;
                    }
                    else { // z axis
                        outfile << intersection.x << " " << intersection.y << " " << intersection.z << " " << 255 <<
                        " " << 0 << " " << 0 << std::endl;
                    }
                }
                outfile.close();
            }
        }
        if (intersection_points.size() > 1) {
            std::cout << numIntersectionsOnEdge << " Interseções na mesma aresta" << std::endl;
            std::ofstream anomalyfile;
            anomalyfile.open("../../anomaly_intersection.ply", std::ios::app);
            for (auto point = intersection_points.begin(); point != intersection_points.end() ; ++point) {
                if (i < 4) { // x axis
                    anomalyfile << point->x << " " << point->y << " " << point->z << " " << 255 << " " <<
                    255 << " " << 0 << std::endl;
                }
                else {
                    if (i < 8) { // y axis
                        anomalyfile << point->x << " " << point->y << " " << point->z << " " << 0 <<
                        " " << 255 << " " << 255 << std::endl;
                    }
                    else { // z axis
                        anomalyfile << point->x << " " << point->y << " " << point->z << " " << 255 <<
                        " " << 0 << " " << 0 << std::endl;
                    }
                }
            }
            anomalyfile.close();
            /*if (numIntersectionsOnEdge%2 == 0){
                vecsigns[c1] = MATERIAL_UNKNOWN;
                vecsigns[c2] = MATERIAL_UNKNOWN;
            }*/
            int mindistance = 9999999;
            int nearindex = -1;
            for (int j = 0; j < intersection_points.size(); ++j) {
                int distance = glm::distance(p1, intersection_points[j]);
                if ( distance < mindistance){
                    nearindex = j;
                    mindistance = distance;
                }
            }
            vecsigns[c1] = glm::dot((p1 - intersection_points[nearindex]), normals[nearindex]) < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
            mindistance = 9999999;
            nearindex = -1;
            for (int j = 0; j < intersection_points.size(); ++j) {
                int distance = glm::distance(p2, intersection_points[j]);
                if ( distance < mindistance){
                    nearindex = j;
                    mindistance = distance;
                }
            }
            vecsigns[c2] = glm::dot((p2 - intersection_points[nearindex]), normals[nearindex]) < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
            /*if (intersection_points.size()%2 == 0) {
                vecsigns[c2] = vecsigns[c1];
            }
            else {
                vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
            }*/

            /*if (numIntersectionsOnEdge == 2){
                //vec3 middle = (intersection_points[0] + intersection_points[1])/2;
                int nearmost_index = glm::distance(p1, intersection_points[0]) < glm::distance(p1, intersection_points[1]) ? 0 : 1;

                vecsigns[c2] = vecsigns[c1];
            }*/
        }
        if (vecsigns[c1] != vecsigns[c2]){
            for (int j = 0; j < intersection_points.size(); ++j) {
                vec3 n = normals[j];
                vec3 v = intersection_points[j];
                qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
                averageNormal += n;
                hasIntersection = true;
                edgedata.hasIntersection = true;
                edgedata.intersection = intersection;
                edgedata.normal = n;
            }
        }

        OctreeNode::edgepool[edgehash] = edgedata;
        std::string vertexhash = hashvertex(p1);
        if (OctreeNode::vertexpool.count(vertexhash) == 0){
            OctreeNode::vertexpool[vertexhash] = vecsigns[c1];
        }
        else {
            if (OctreeNode::vertexpool[vertexhash] == MATERIAL_UNKNOWN /*|| vecsigns[c2] == MATERIAL_SOLID*/){
                OctreeNode::vertexpool[vertexhash] = vecsigns[c1];
            }
        }
        vertexhash = hashvertex(p2);
        if (OctreeNode::vertexpool.count(vertexhash) == 0){
            OctreeNode::vertexpool[vertexhash] = vecsigns[c2];
        }
        else {
            if (OctreeNode::vertexpool[vertexhash] == MATERIAL_UNKNOWN /*|| vecsigns[c2] == MATERIAL_SOLID*/){
                OctreeNode::vertexpool[vertexhash] = vecsigns[c2];
            }
        }

    }
    if (!hasIntersection){
        // voxel is full inside or outside the volume
        delete leaf;
        return nullptr;
    }

    //std::cout << std::endl;
    int ground_corners = 0;
    bool checksigns = true;
    while(checksigns) {
        checksigns = false;
        for (size_t i = 0; i < 8; i++) {
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



    std::ofstream gridfile;
    //gridfile.open("/home/hallpaz/Workspace/dual_contouring_experiments/grid_color_updated.ply", std::ios::app);
    gridfile.open("../../grid_color_updated.ply", std::ios::app);
    for (size_t i = 0; i < 8; i++) {
        corners |= (vecsigns[i] << i);

        const vec3 cornerPos = leaf->min + CHILD_MIN_OFFSETS[i]*leaf->size;
//        const float density = Density_Func(vec3(cornerPos));
//        const int material = density < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
//        ground_corners |= (material << i);
        if (vecsigns[i] == MATERIAL_SOLID) {
            gridfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 128 << " " << 255 << std::endl;
        }

        if (vecsigns[i] == MATERIAL_AIR) {
            gridfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 128 << " " << 128 << " " << 64 << std::endl;
        }

        if (vecsigns[i] == MATERIAL_UNKNOWN) {
            gridfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
        }
//        const float density = Density_Func(vec3(cornerPos));
//        const int material = density < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
//        ground_corners |= (material << i);
    }
    gridfile.close();
//    if (corners != ground_corners) {
//        std::cout << " edgecount: " << edgecount << " corners: " << corners << " ground: " << ground_corners << std::endl;
//        for (int i = 0; i < 8; ++i) {
//            std::cout << ((corners >> i) & 1) << " x " << (vecsigns[i]) << " x " << ((ground_corners >> i) & 1) << std::endl;
//        }
//        //exit(1);
//    }

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

    drawInfo->averageNormal = glm::normalize(averageNormal / (float)edgecount);
    drawInfo->corners = corners;

    leaf->type = Node_Leaf;
    leaf->drawInfo = drawInfo;

    return leaf;
}

RelativePosition OctreeNode::vertexRelativePosition(const Vertex &vertex) {
    /*We ignore the case when the vertex is exactly on the cell edge (we consider it outside) */

    vec3 maxvertex = this->min + (CHILD_MIN_OFFSETS[7] * this->size);
//    std::cout <<"min: (" << min.x << ", " << min.y << ", " << min.z << ")" << std::endl;
//    std::cout <<"max: (" << maxvertex.x << ", " << maxvertex.y << ", " << maxvertex.z << ")" << std::endl;
//    std::cout <<"vertex: (" << vertex.position.x << ", " << vertex.position.y << ", " << vertex.position.z << ")" << std::endl;
    if ((vertex.position.x > this->min.x) && (vertex.position.x < maxvertex.x)
            && (vertex.position.y > this->min.y) && (vertex.position.y < maxvertex.y)
            && (vertex.position.z > this->min.z) && (vertex.position.z < maxvertex.z)) {
        return INSIDE;
    }
    return OUTSIDE;
}

RelativePosition OctreeNode::triangleRelativePosition(const Vertex &a, const Vertex &b, const Vertex &c) {
    vec3 max = min + size*CHILD_MIN_OFFSETS[7];
    if (((a.position.x > max.x) && (b.position.x > max.x) && (c.position.x > max.x))
        || ((a.position.x < min.x) && (b.position.x < min.x) && (c.position.x < min.x))){
        return OUTSIDE;
    }
    if (((a.position.y > max.y) && (b.position.y > max.y) && (c.position.y > max.y))
        || ((a.position.y < min.y) && (b.position.y < min.y) && (c.position.y < min.y))){
        return OUTSIDE;
    }
    if (((a.position.z > max.z) && (b.position.z > max.z) && (c.position.z > max.z))
        || ((a.position.z < min.z) && (b.position.z < min.z) && (c.position.z < min.z))){
        return OUTSIDE;
    }

    int numVerticesInside = 0;
    if (this->vertexRelativePosition(a) == INSIDE){
        ++numVerticesInside;
    }
    if (this->vertexRelativePosition(b) == INSIDE){
        ++numVerticesInside;
    }
    if (this->vertexRelativePosition(c) == INSIDE){
        ++numVerticesInside;
    }

    if (numVerticesInside == 3) return INSIDE;


    return CROSSING;

}




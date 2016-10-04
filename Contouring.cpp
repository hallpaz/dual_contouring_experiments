//
// Created by Hallison da Paz on 02/09/2016.
//

#include "Contouring.h"
#include "Constants.h"

#include "glm/glm.hpp"

using glm::vec3;

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
        svd::Vec3 qefPosition;
        svd::QefSolver qef;
        qef.setData(d->qef);

        qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
        d->position = vec3(qefPosition.x, qefPosition.y, qefPosition.z);
        d->qef = qef.getData();

        const vec3 min = vec3(node->min);
        const vec3 max = vec3(node->min + vec3(node->size));
        if (d->position.x < min.x || d->position.x > max.x ||
            d->position.y < min.y || d->position.y > max.y ||
            d->position.z < min.z || d->position.z > max.z)
        {
            const auto& mp = qef.getMassPoint();
            d->position = vec3(mp.x, mp.y, mp.z);
        }
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
            indexBuffer.push_back(Triangle{indices[0], indices[1], indices[3]});
            indexBuffer.push_back(Triangle{indices[0], indices[3], indices[2]});
        }
        else
        {
            indexBuffer.push_back(Triangle{indices[0], indices[3], indices[1]});
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
                if (node[j]->type != Node_Internal/*node[j]->type == Node_Leaf || node[j]->type == Node_Pseudo*/)
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
                if (node[order[j]]->type != Node_Internal/*node[order[j]]->type == Node_Leaf || node[order[j]]->type == Node_Pseudo*/)
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
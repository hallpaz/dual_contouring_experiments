//
// Created by Hallison da Paz on 02/09/2016.
//

#include <fstream>
#include "Contouring.h"
#include "Constants.h"

#include "glm/glm.hpp"

using glm::vec3;


// Declarations
void GenerateVertexIndices(OctreeNode* node, VertexBuffer& vertexBuffer);
void ContourProcessEdge(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer);
void ContourEdgeProc(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer);
void ContourFaceProc(OctreeNode* node[2], int dir, IndexBuffer& indexBuffer);
void ContourCellProc(OctreeNode* node, IndexBuffer& indexBuffer);

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

    if (node->type != NODE_LEAF)
    {
        for (int i = 0; i < NUM_CHILDREN; i++)
        {
            GenerateVertexIndices(node->children[i], vertexBuffer);
        }
    }

    if (node->type != NODE_INTERNAL)
    {
        OctreeDrawInfo* d = node->drawInfo;
        if (!d)
        {
            std::cout << "Error! Could not add vertex! DrawInfo is not available" << std::endl;
            exit(EXIT_FAILURE);
        }

        Vertex vertex;
        d->index = vertexBuffer.size();
        d->averageNormal = glm::normalize(d->averageNormal);
        svd::Vec3 qefPosition;
        svd::QefSolver qef;
        qef.setData(d->qef);
        qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
        d->position = vec3(qefPosition.x, qefPosition.y, qefPosition.z);

        const vec3 min = vec3(node->min);
        const vec3 max = vec3(node->min + vec3(node->size));
        if (d->position.x < min.x || d->position.x > max.x ||
            d->position.y < min.y || d->position.y > max.y ||
            d->position.z < min.z || d->position.z > max.z)
        {
            const auto& mp = qef.getMassPoint();
            d->position = vec3(mp.x, mp.y, mp.z);
            vertex.color = glm::uvec3(255, 0, 0);
            if (node->is_border){
                vertex.color = glm::uvec3(0, 255, 0);
            }
        }
        else
        {
            vertex.color = glm::uvec3(128, 128, 128);
            if (node->is_border){
                vertex.color = glm::uvec3(0, 0, 255);
            }
        }
        vertex.position = d->position;
        //EVERY POINT IN THE MIDDLE OF THE CELL
        //d->position = vec3(node->min + vec3(node->size/2));
        //d->position = glm::normalize(vec3(node->min + vec3(node->size/2)))*4.0f;
        // ----------------------------------------------
        vertexBuffer.push_back(/*Vertex(d->position)*/vertex);

        // DEBUG --------------------------------------------------------------------
        std::ofstream interiorfile, exteriorfile;
        interiorfile.open("../subproducts/interior_color_updated.ply", std::ios::app);
        exteriorfile.open("../subproducts/exterior_color_updated.ply", std::ios::app);
        for (int i = 0; i < NUM_CHILDREN; ++i) {
            const vec3 cornerPos = node->get_vertex(i);
            int sign = (node->drawInfo->corners >> i) & 1;
            if (sign == MATERIAL_SOLID) {
                interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 128 << " " << 255 << std::endl;
            }

            if (sign == MATERIAL_AIR) {
                exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 64 << " " << 255 << " " << 64 << std::endl;
            }

            if (sign == MATERIAL_UNKNOWN) {
                interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
                exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
            }

        }
        interiorfile.close();
        exteriorfile.close();
        // DEBUG --------------------------------------------------------------------//
    }
}

// ----------------------------------------------------------------------------

void ContourProcessEdge(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer)
{
    Real minSize = 1000000.f;		// arbitrary big number
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
        if (indices[0] != indices[3])
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
}

// ----------------------------------------------------------------------------

void ContourEdgeProc(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer)
{
    if (!node[0] || !node[1] || !node[2] || !node[3])
    {
        return;
    }

    if (node[0]->type != NODE_INTERNAL &&
        node[1]->type != NODE_INTERNAL &&
        node[2]->type != NODE_INTERNAL &&
        node[3]->type != NODE_INTERNAL)
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
                if (node[j]->type != NODE_INTERNAL/*node[j]->type == NODE_LEAF || node[j]->type == NODE_PSEUDO*/)
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

    if (node[0]->type == NODE_INTERNAL ||
        node[1]->type == NODE_INTERNAL)
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
                if (node[j]->type != NODE_INTERNAL)
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
                if (node[order[j]]->type != NODE_INTERNAL/*node[order[j]]->type == NODE_LEAF || node[order[j]]->type == NODE_PSEUDO*/)
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
    if (node == nullptr)
    {
        return;
    }

    if (node->type == NODE_INTERNAL)
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
//
// Created by Hallison da Paz on 02/09/2016.
//

#ifndef DUAL_CONTOURING_EXPERIMENTS_CONTOURING_H
#define DUAL_CONTOURING_EXPERIMENTS_CONTOURING_H

#include "octree.h"


void GenerateMeshFromOctree(OctreeNode* node, VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer);
void GenerateVertexIndices(OctreeNode* node, VertexBuffer& vertexBuffer);
void ContourProcessEdge(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer);
void ContourEdgeProc(OctreeNode* node[4], int dir, IndexBuffer& indexBuffer);
void ContourFaceProc(OctreeNode* node[2], int dir, IndexBuffer& indexBuffer);
void ContourCellProc(OctreeNode* node, IndexBuffer& indexBuffer);


#endif //DUAL_CONTOURING_EXPERIMENTS_CONTOURING_H

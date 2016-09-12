//
// Created by Hallison da Paz on 11/09/2016.
//

#include "Tests.h"
#include "Constants.h"
#include "Utils.h"
#include "glm/glm.hpp"

using glm::vec3;
using std::string;

bool Tests::validate_cells_signs(OctreeNode *node, std::unordered_map<std::string, int> &vertexpool, int &num_incorrect)
{
    if (node == nullptr)
    {
        return true;
    }
    if (node->type == Node_Internal)
    {
        for (int i = 0; i < 8; ++i) {
            validate_cells_signs(node->children[i], vertexpool, num_incorrect);
        }
    }
    else
    {
        int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                           MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};

        for (int i = 0; i < 8; ++i) {
            vec3 vertex = vec3(node->min + node->size*CHILD_MIN_OFFSETS[i]);
            string vertexhash = hashvertex(vertex);
            if (vertexpool.count(vertexhash) == 0)
            {
                std::cout << "Vertex NOT REGISTERED" << std::endl;
                ++num_incorrect;
                return false;
            }
            vecsigns[i] = vertexpool[vertexhash];
            if (vecsigns[i] == -1)
            {
                std::cout << "Vertex Sign Inconsistency MATERIAL UNKNOWN" << std::endl;
                ++num_incorrect;
                return false;
            }
        }
        int local_corners = 0;
        for (size_t i = 0; i < 8; i++)
        {   //encode the signs to the corners variable to save memory
            local_corners |= (vecsigns[i] << i);
        }
        if (local_corners != node->drawInfo->corners){
            std::cout << local_corners << " | " << node->drawInfo->corners << std::endl;
            ++num_incorrect;
        }
    }
    if (num_incorrect == 0)
    {
        return true;
    }
    return false;
}


bool Tests::validate_vertices_map(std::unordered_map<std::string, int> &vertexpool)
{
    int invalids = 0;
    for (auto it = vertexpool.begin(); it != vertexpool.end(); ++it)
    {
        if (it->second == MATERIAL_UNKNOWN)
        {
            invalids++;
        }
    }
    std::cout << "Total Signs: " << vertexpool.size() << "| Number of Invalid Signs: " << invalids << std::endl;
    if (invalids == 0) return true;
    return false;
}
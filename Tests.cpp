//
// Created by Hallison da Paz on 11/09/2016.
//

#include "Tests.h"
#include "Constants.h"
#include "Utils.h"
#include "glm/glm.hpp"
#include "NormalsEstimator.h"
#include "Reconstruction.h"
#include "Contouring.h"

using glm::vec3;
using std::string;

bool Tests::validate_cells_signs(OctreeNode *node, std::unordered_map<std::string, int> &vertexpool, int &num_incorrect)
{
    if (node == nullptr)
    {
        return true;
    }
    if (node->type == NODE_INTERNAL)
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


string folder_for_parameters(string basename, int pieces, bool should_simplify)
{
    std::stringstream folderpath;
    folderpath << "../experiments/" << basename + "/";
    if (should_simplify)
    {
        folderpath << "simplified/";
    }
    folderpath << pieces + "/";
    return folderpath.str();
}

void Tests::reconstruct_pieces(string input_folder, string basename, int num_pieces, int height, bool should_simplify){
    OctreeNode* root = nullptr;
    const string EXTENSION = ".off";

    string folder_name = folder_for_parameters(basename, num_pieces, should_simplify);

    string bboxfilename = input_folder + basename + EXTENSION;
    DefaultMesh myMesh;
    OpenMesh::IO::read_mesh(myMesh, bboxfilename);
    // compute bounding box
    DefaultMesh::Point bb_min, bb_max;
    auto v_it = myMesh.vertices_begin();
    bb_min = bb_max = myMesh.point(*v_it);
    float big = -1;
    float small = 10000;
    for (; v_it != myMesh.vertices_end(); ++v_it)
    {
        bb_min.minimize(myMesh.point(*v_it));
        bb_max.maximize(myMesh.point(*v_it));

        float value =  myMesh.point(*v_it).length();
        if (value > big){
            big = value;
        }
        if (value < small){
            small = value;
        }
    }
    std::cout << big << " " << small << std::endl;
    float octreeSize = (bb_max - bb_min).max();
    std::cout << "Min: (" << bb_min[0] << ", " << bb_min[1] << ", " << bb_min[2] << ") " << "Size: " << octreeSize << std::endl;

    myMesh.request_vertex_status();
    myMesh.request_edge_status();
    myMesh.request_face_status();
    NormalsEstimator::compute_better_normals(myMesh);

    float simpthreshold = octreeSize/1000.0;

    VertexBuffer vertices;
    IndexBuffer indices;

    std::cout << "TEST: will start build" << std::endl;
    std::vector<string> filenames;
    for (int i = 0; i < num_pieces; ++i) {
        filenames.push_back(basename + i + EXTENSION);
    }
    root = Fusion::octree_from_samples(glm::vec3(bb_min[0], bb_min[1], bb_min[2]) - vec3(0.1), octreeSize+0.2/**1.1*/, height, filenames);
    if (should_simplify){
        root = SimplifyOctree(root, simpthreshold);
    }

    std::cout << "TEST: will start mesh generation" << std::endl;
    GenerateMeshFromOctree(root, vertices, indices);
    std::cout << vertices.size() << std::endl;
    std::cout << indices.size() << std::endl;
    std::stringstream filepath;
    filepath << folder_name << basename << height << EXTENSION;
    write_OFF(filepath.str(), vertices, indices);
    std::cout << "Generated mesh" << std::endl;


    return;
}
#include <iostream>
#include <sstream>
#include <vector>
#include "Utils.h"
#include "octree.h"

#include "NormalsEstimator.h"
#include "Reconstruction.h"
#include "Contouring.h"


using namespace std;
using glm::vec3;


// ----------------------------------------------------------------------------

int main(int argc, char** argv)
{
    OctreeNode* root = nullptr;
    const int height = 9;

    string folder_name = "../";
    //std::cout << "Input File Name" << endl;
    string inputfilename, outputfilename;
    //std::cin >> inputfilename;
    std::cout <<"Output File Name" << endl;
    std::cin >> outputfilename;

    DefaultMesh myMesh;
//    OpenMesh::IO::read_mesh(myMesh, folder_name + inputfilename);
    OpenMesh::IO::read_mesh(myMesh, "../models/sphere8.off");
    // compute bounding box
    DefaultMesh::Point bb_min, bb_max;
    auto v_it = myMesh.vertices_begin();
    bb_min = bb_max = myMesh.point(*v_it);
    for (; v_it != myMesh.vertices_end(); ++v_it)
    {
        bb_min.minimize(myMesh.point(*v_it));
        bb_max.maximize(myMesh.point(*v_it));
    }
    float octreeSize = (bb_max - bb_min).max();
    std::cout << "Min: (" << bb_min[0] << ", " << bb_min[1] << ", " << bb_min[2] << ") " << "Size: " << octreeSize << std::endl;

    myMesh.request_vertex_status();
    myMesh.request_edge_status();
    myMesh.request_face_status();
    NormalsEstimator::compute_better_normals(myMesh);

//    float simpthreshold = octreeSize/10.0;

    cout << "Generating mesh with octreeSize: " << octreeSize << "\n" << endl;

    VertexBuffer vertices;
    IndexBuffer indices;

    cout << "MAIN: will start build" << endl;
    //root = BuildOctree(glm::vec3(-octreeSize / 2), octreeSize, height, simpthreshold);
//    root = BuildOctreeFromMesh(minPoint, octreeSize, height, simpthreshold, testVertices, testIndices);
//    root = BuildOctreeFromOpenMesh(glm::vec3(bb_min[0], bb_min[1], bb_min[2]) - vec3(0.1), octreeSize*1.1, height, myMesh);

    std::vector<string> filenames { "../models/divided/esfera3.off", "../models/divided/esfera0.off"};//, "../models/divided/esfera2.off", "../models/divided/esfera1.off"};
    //std::vector<string> filenames {"../models/sphere8.off"};
    root = Fusion::octree_from_samples(glm::vec3(bb_min[0], bb_min[1], bb_min[2]) - vec3(0.1), octreeSize*1.1, height, filenames);

    /*if (!Tests::validate_vertices_map(OctreeNode::vertexpool))
    {
        exit(99);
    }
    int num_incorrect = 0;
    if (Tests::validate_cells_signs(root, OctreeNode::vertexpool, num_incorrect))
    {
        cout << "ALL Signs are correct" << endl;
    }
    else
    {
        cout << "There are " << num_incorrect << " signs incorrect" << endl;
    }*/

//    std::cout << "Octree.cpp: will simplify nodes" << std::endl;
//    root = SimplifyOctree(root, simpthreshold);
//    std::cout << "Octree.cpp: did simplify nodes" << std::endl;
    cout << "MAIN: will start mesh generation" << endl;
    GenerateMeshFromOctree(root, vertices, indices);
    cout << vertices.size() << endl;
    cout << indices.size() << endl;
    std::stringstream filepath;
    filepath << folder_name << outputfilename << height << ".off";
    write_OFF(filepath.str(), vertices, indices);
    printf("Generated mesh\n\n");


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

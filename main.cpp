#include <iostream>
#include <sstream>
#include <vector>
#include "Utils.h"
#include "octree.h"

#include "NormalsEstimator.h"


using namespace std;
using glm::vec3;
// ----------------------------------------------------------------------------

int main(int argc, char** argv)
{
    bool IMPA = false;
    OctreeNode* root = nullptr;
    const int height = 6;


    string folder_name = "../";
    std::cout << "Input File Name" << endl;
    string inputfilename, outputfilename;
    std::cin >> inputfilename;
    std::cout <<"Output File Name" << endl;
    std::cin >> outputfilename;


    DefaultMesh myMesh;
    OpenMesh::IO::read_mesh(myMesh, folder_name + inputfilename);
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
//    OctreeNode *root = BuildOctreeFromOpenMesh(glm::vec3(bb_min[0], bb_min[1], bb_min[2]), octreeSize, height, simpthreshold, myMesh);

    float simpthreshold = octreeSize/100000.0;

    cout << "Generating mesh with octreeSize: " << octreeSize << "\n" << endl;

    VertexBuffer vertices;
    IndexBuffer indices;

    cout << "MAIN: will start build" << endl;
    //root = BuildOctree(glm::vec3(-octreeSize / 2), octreeSize, height, simpthreshold);
//    root = BuildOctreeFromMesh(minPoint, octreeSize, height, simpthreshold, testVertices, testIndices);
    root = BuildOctreeFromOpenMesh(glm::vec3(bb_min[0], bb_min[1], bb_min[2]) - vec3(0.1), octreeSize*1.1, height, myMesh);
    std::cout << "Octree.cpp: will simplify nodes" << std::endl;
    /*root = SimplifyOctree(root, simpthreshold);
    std::cout << "Octree.cpp: did simplify nodes" << std::endl;
    cout << "MAIN: will start mesh generation" << endl;*/
    GenerateMeshFromOctree(root, vertices, indices);
    cout << vertices.size() << endl;
    cout << indices.size() << endl;
    std::stringstream filepath;
    filepath << folder_name << outputfilename << height << ".off";
    write_OFF(filepath.str(), vertices, indices);
    printf("Generated mesh\n\n");
    simpthreshold = simpthreshold/10.0;


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

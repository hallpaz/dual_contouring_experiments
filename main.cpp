#include <iostream>
#include <sstream>
#include "Utils.h"
#include "glm/glm.hpp"
#include "octree.h"


using namespace std;
// ----------------------------------------------------------------------------




int main(int argc, char** argv)
{
    bool IMPA = false;
    OctreeNode* root = nullptr;
    // octreeSize must be a power of two!
    float octreeSize = 4;
    const int height = 6;

    VertexBuffer testVertices;
    IndexBuffer testIndices;

    string folder_name = IMPA ?  "/home/hallpaz/Workspace/dual_contouring_experiments/" : "/Users/hallpaz/Workspace/research/dual_contouring_experiments/";

    octreeSize = read_OFF(testVertices, testIndices, folder_name + "models/cow.off");
    //write_OFF(testVertices, testIndices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/testCube.off");
    //octreeSize = read_OFF(testVertices, testIndices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/debug/sphere51.off");
    std::cout << "num of vertices: " << testVertices.size() << " num of indices: " << testIndices.size() << std::endl;

    float simpthreshold = 0.1;
    for (int i = 1; i < 2; ++i)
    {

        //thresholdIndex = (thresholdIndex + 1) % MAX_THRESHOLDS;
        cout << "Generating mesh with octreeSize: " << octreeSize << "\n" << endl;

        VertexBuffer vertices;
        IndexBuffer indices;

        cout << "MAIN: will start build" << endl;
        //root = BuildOctree(glm::vec3(-octreeSize / 2), octreeSize, height, simpthreshold);
        root = BuildOctreeFromMesh(glm::vec3(-octreeSize / 2), octreeSize, height, simpthreshold, testVertices, testIndices);
        cout << "MAIN: will start mesh generation" << endl;
        GenerateMeshFromOctree(root, vertices, indices);
        cout << vertices.size() << endl;
        cout << indices.size() << endl;
        //write_OFF(vertices, indices, folder_name + "dc6_vase.off");
        std::stringstream filepath;
        filepath << folder_name << "check/simpcowdc" << height << i << octreeSize << ".off";
        write_OFF(vertices, indices, filepath.str());
        printf("Generated mesh\n\n");
        simpthreshold = simpthreshold/10.0;
    }


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

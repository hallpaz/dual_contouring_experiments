#include <iostream>
#include <sstream>
#include "Utils.h"
#include "octree.h"


using namespace std;
// ----------------------------------------------------------------------------




int main(int argc, char** argv)
{
    bool IMPA = false;
    OctreeNode* root = nullptr;
    // octreeSize must be a power of two!
    float octreeSize = 4;
    const int height = 7;

    VertexBuffer testVertices;
    IndexBuffer testIndices;

    string folder_name = IMPA ?  "/home/hallpaz/Workspace/dual_contouring_experiments/" : "/Users/hallpaz/Workspace/research/dual_contouring_experiments/";

    octreeSize = read_OFF(testVertices, testIndices, folder_name + "models/sofa_better.off");
    std::cout << "num of vertices: " << testVertices.size() << " num of indices: " << testIndices.size() << std::endl;

    float simpthreshold = octreeSize/100000.0;
    for (int i = 1; i < 2; ++i)
    {
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
        std::stringstream filepath;
        filepath << folder_name << "check/sofabetterdc_simp" << height << i << ".off";
        write_OFF(vertices, indices, filepath.str());
        printf("Generated mesh\n\n");
        simpthreshold = simpthreshold/10.0;
    }


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

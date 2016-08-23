#include <iostream>
#include <sstream>
#include <vector>
#include "Utils.h"
#include "octree.h"


using namespace std;
// ----------------------------------------------------------------------------

void computePerVertexNormals(VertexBuffer &vertices, IndexBuffer &indices, VertexBuffer &normals)
{
    std::vector<std::vector<int>> inclusionList(vertices.size());
    for (int i = 0; i < indices.size(); ++i) {

        inclusionList[indices[i].a].push_back(i);
        inclusionList[indices[i].b].push_back(i);
        inclusionList[indices[i].c].push_back(i);
    }

    for (int j = 0; j < inclusionList.size(); ++j) {
        for (std::vector<int>::iterator triangle = inclusionList[j].begin(); triangle != inclusionList[j].end() ; ++triangle) {
            
        }
    }
}


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
    std::cout << "Input File Name" << endl;
    string inputfilename, outputfilename;
    std::cin >> inputfilename;
    std::cout <<"Output File Name" << endl;
    std::cin >> outputfilename;
    vec3 minPoint = vec3(0, 0, 0);
    octreeSize = read_OFF(folder_name + inputfilename, testVertices, testIndices, minPoint);
    std::cout << "num of vertices: " << testVertices.size() << " num of indices: " << testIndices.size() << std::endl;

    float simpthreshold = octreeSize/100000.0;
    for (int i = 1; i < 2; ++i)
    {
        cout << "Generating mesh with octreeSize: " << octreeSize << "\n" << endl;

        VertexBuffer vertices;
        IndexBuffer indices;

        cout << "MAIN: will start build" << endl;
        //root = BuildOctree(glm::vec3(-octreeSize / 2), octreeSize, height, simpthreshold);
        //root = BuildOctreeFromMesh(glm::vec3(-octreeSize / 2), octreeSize, height, simpthreshold, testVertices, testIndices);
        root = BuildOctreeFromMesh(minPoint, octreeSize, height, simpthreshold, testVertices, testIndices);
        cout << "MAIN: will start mesh generation" << endl;
        GenerateMeshFromOctree(root, vertices, indices);
        cout << vertices.size() << endl;
        cout << indices.size() << endl;
        std::stringstream filepath;
        filepath << folder_name << outputfilename << height << i << ".off";
        write_OFF(filepath.str(), vertices, indices);
        printf("Generated mesh\n\n");
        simpthreshold = simpthreshold/10.0;
    }


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

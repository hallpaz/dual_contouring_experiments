#include <iostream>

#include "Utils.h"
#include "glm/glm.hpp"
#include "octree.h"

using namespace std;
// ----------------------------------------------------------------------------

/*void write_OFF(string filename, VertexBuffer &vertices, IndexBuffer &indices)
{
    cout << "Entrou no write" << endl;
    std::ofstream meshfile("ninja.txt");
    cout << "WTF????????" << endl;
    //meshfile.open(filename.c_str());
    cout << "pqssou" << endl;
    if(meshfile.is_open()) {
        std::cout << "Arquivo aberto: " << filename << std::endl;

        meshfile << "OFF" << std::endl;
        meshfile << vertices.size() << " " << indices.size() << " 0" << std::endl;

        *//*for ( std::vector<glm::vec3>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) {
            meshfile << v_it->x << " " << v_it->y << " " << v_it->z << std::endl;
        }

        for ( std::vector<int>::iterator f_it = indices.begin(); f_it != indices.end(); f_it+=3) {
            meshfile << "3" << " " << *f_it << " " << *(f_it+1) << " " << *(f_it+2) << std::endl;
        }*//*
        meshfile.close();
    }
    else {
        cout << "Unable to open the file" << endl;
    }

}*/

int main(int argc, char** argv)
{

    OctreeNode* root = nullptr;
    // octreeSize must be a power of two!
    const int octreeSize = 3.0;
    const int height = 4;

    bool refreshMesh = true;
    int thresholdIndex = 0;

    if (refreshMesh)
    {
        refreshMesh = false;
        //thresholdIndex = (thresholdIndex + 1) % MAX_THRESHOLDS;
        cout << "Generating mesh with error threshold...\n" << endl;

        VertexBuffer vertices;
        IndexBuffer indices;

        cout << "MAIN: will start build" << endl;
        root = BuildOctree(glm::vec3(-octreeSize / 2), octreeSize, height, 0.0000000001);
        cout << "MAIN: will start mesh generation" << endl;
        GenerateMeshFromOctree(root, vertices, indices);
        cout << vertices.size() << endl;
        cout << indices.size() << endl;
        write_OFF(vertices, indices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/torusH.off");
        printf("Generated mesh\n\n");
    }


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

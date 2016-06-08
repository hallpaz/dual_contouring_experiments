#include <iostream>
#include <sstream>
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

void createMesh(VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer){
    int paralelos = 20;
    int meridianos = 40;

    float theta = 0.005;
    float phi = 0.00;

    int j = 0; int i = 0;
    float u, v, depth;

    for(theta = 0.00005; theta < M_PI; theta += (M_PI - 0.0001)/(paralelos-1)){
        for(phi = 0.0; phi < 2*M_PI; phi += 2*M_PI/meridianos){
            //Position coordinates
            depth = 16.0f;

            u = phi/(2*M_PI);
            v = 1.0 - theta/M_PI;

            vertexBuffer.push_back(Vertex(glm::vec3(depth*sin(theta)*sin(phi), depth*cos(theta), depth*sin(theta)*cos(phi))));
            ++i;
        }
    }
    //compute triangulation
    double distanceThreshold = 0.2;
    double mindist = 100.0, maxdist = -1.0;
    for(i = 0; i < paralelos-2; ++i){
        for (j = 0; j < meridianos-1; ++j) {

            indexBuffer.push_back( { i*(meridianos+1) + j, (i+1)*(meridianos+1) +j, i*(meridianos+1) + j+1  });


            indexBuffer.push_back( {(i+1)*(meridianos+1) + j+1, i*(meridianos+1) + j, i*(meridianos+1) + j+1});
        }
    }
}

void createCylinder(VertexBuffer& vertexBuffer, IndexBuffer& indexBuffer){
    int paralelos = 40;
    int meridianos = 40;



    int j = 0; int i = 0;
    float demiheight = 8.0f;

    for(float theta = -demiheight; theta < demiheight; theta += (2*demiheight)/(paralelos)){
        for(float phi = 0.0; phi < 2*M_PI; phi += 2*M_PI/meridianos){
            //Position coordinates
            //height = 16.0f;

            vertexBuffer.push_back(Vertex(glm::vec3(demiheight*cos(phi), demiheight*sin(phi), theta)));
            ++i;
        }
    }
    //compute triangulation
    for(i = 0; i < paralelos-1; ++i){
        for (j = 0; j < meridianos; ++j) {

            indexBuffer.push_back( { i*(meridianos) + j, i*(meridianos) + j+1, (i+1)*(meridianos) +j  });


            indexBuffer.push_back( {(i+1)*(meridianos) + j+1, (i+1)*(meridianos) + j, i*(meridianos) + j+1});
        }
    }
}



int main(int argc, char** argv)
{

    OctreeNode* root = nullptr;
    // octreeSize must be a power of two!
    float octreeSize = 4;
    const int height = 6;

    bool refreshMesh = true;
    int thresholdIndex = 0;

    VertexBuffer testVertices;
    IndexBuffer testIndices;
    //createCylinder(testVertices, testIndices);

    octreeSize = read_OFF(testVertices, testIndices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/models/vase_mesh_simp.off");
    //write_OFF(testVertices, testIndices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/testCube.off");
    //octreeSize = read_OFF(testVertices, testIndices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/debug/sphere51.off");
    std::cout << "num of vertices: " << testVertices.size() << " num of indices: " << testIndices.size() << std::endl;

    /*if (refreshMesh)
    {
        refreshMesh = false;
        //thresholdIndex = (thresholdIndex + 1) % MAX_THRESHOLDS;
        cout << "Generating mesh with octreeSize: " << octreeSize << "\n" << endl;

        VertexBuffer vertices;
        IndexBuffer indices;

        cout << "MAIN: will start build" << endl;
        root = BuildOctree(glm::vec3(-octreeSize / 2), octreeSize, height, 0.0001);
        //root = BuildOctreeFromMesh(glm::vec3(-octreeSize / 2), octreeSize, height, 0.001, testVertices, testIndices);
        cout << "MAIN: will start mesh generation" << endl;
        GenerateMeshFromOctree(root, vertices, indices);
        cout << vertices.size() << endl;
        cout << indices.size() << endl;
        //write_OFF(vertices, indices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/dc6_teddy1596.off");
        write_OFF(vertices, indices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/sphere6.off");
        printf("Generated mesh\n\n");
    }*/

    float simpthreshold = 0.001;
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
        write_OFF(vertices, indices, "/Users/hallpaz/Workspace/research/dual_contouring_experiments/dc6_vase.off");
        /*std::stringstream filepath;
        filepath << "/Users/hallpaz/Workspace/research/dual_contouring_experiments/sphere_poly" << height << i << octreeSize << ".off";
        write_OFF(vertices, indices, filepath.str());*/
        printf("Generated mesh\n\n");
        simpthreshold = simpthreshold/10.0;
    }


    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

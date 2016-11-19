#include <iostream>
#include <sstream>
#include <vector>
#include "Utils.h"
#include "old/octree.h"

#include "old/NormalsEstimator.h"
#include "old/Reconstruction.h"
#include "Contouring.h"
#include "Octree.h"


using namespace std;
using glm::vec3;


// ----------------------------------------------------------------------------

Real compute_boundingbox(DefaultMesh &mesh, DefaultMesh::Point &bb_min);

int main(int argc, char** argv)
{
    const int height = 9;

    string folder_name = "../";
    string inputfilename, outputfilename;
    std::cout <<"Output File Name" << endl;
    std::cin >> outputfilename;

    DefaultMesh myMesh;
    OpenMesh::IO::read_mesh(myMesh, "../models/analytic/sphere_lowpoly.off");
    // compute bounding box
    DefaultMesh::Point bb_min;
    Real octreeSize = compute_boundingbox(myMesh, bb_min);

    Octree sphere_octree(openmesh_to_glm(bb_min) - vec3(0.1), octreeSize*1.1, height, myMesh);

    VertexBuffer vertices;
    IndexBuffer indices;
    GenerateMeshFromOctree(sphere_octree.root, vertices, indices);

    std::stringstream filepath;
    filepath << folder_name << outputfilename << height << ".off";
    write_OFF(filepath.str(), vertices, indices);
    return EXIT_SUCCESS;
}

// ----------------------------------------------------------------------------

Real compute_boundingbox(DefaultMesh &mesh, DefaultMesh::Point &bb_min)
{
    auto v_it = mesh.vertices_begin();
    DefaultMesh::Point bb_max = mesh.point(*v_it);
    bb_min = bb_max;
    for (; v_it != mesh.vertices_end(); ++v_it)
    {
        bb_min.minimize(mesh.point(*v_it));
        bb_max.maximize(mesh.point(*v_it));
    }
    Real octreeSize = (bb_max - bb_min).max();
    std::cout << "Min: (" << bb_min[0] << ", " << bb_min[1] << ", " << bb_min[2] << ") " << "Size: " << octreeSize << std::endl;
    return octreeSize;
}
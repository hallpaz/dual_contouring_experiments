#include <iostream>
#include <sstream>
#include <vector>
#include "Utils.h"
#include "old/octree.h"

#include "old/NormalsEstimator.h"
#include "Reconstruction.h"
#include "Contouring.h"
#include "Octree.h"


using namespace std;
using glm::vec3;


// ----------------------------------------------------------------------------

Real compute_boundingbox(DefaultMesh &mesh, DefaultMesh::Point &bb_min);

int main(int argc, char** argv)
{
    const int height = 15;
    int dist = 16;

    string folder_name = "../";
    string inputfilename, outputfilename;
    std::cout <<"Output File Name" << endl;
    std::cin >> outputfilename;

    DefaultMesh myMesh;
    //OpenMesh::IO::read_mesh(myMesh, "../models/analytic/sphere_lowpoly.off");
    OpenMesh::IO::read_mesh(myMesh, "../models/divided/vase_antonina.off");
    // compute bounding box
    DefaultMesh::Point bb_min;
    Real octreeSize = compute_boundingbox(myMesh, bb_min);

    //dist = (int) octreeSize;

    //OpenMesh::IO::read_mesh(myMesh, "../models/analytic/sphere_low2.off");
    //Octree sphere_octree(openmesh_to_glm(bb_min) - vec3(0.1), octreeSize*1.1, height, myMesh);
    /*std::vector<std::string> filenames = {"../models/analytic/sphere_lowPZ.off", //"../models/analytic/sphere_lowMZ.off",
                                          "../models/analytic/sphere_lowPX.off", "../models/analytic/sphere_lowMX.off",
                                          "../models/analytic/sphere_lowPY.off", "../models/analytic/sphere_lowMY.off"};*/

    std::vector<std::string> filenames = {"../models/divided/vase_antoninaPZ.off",
                                "../models/divided/vase_antoninaPX.off", "../models/divided/vase_antoninaMX.off",
                                "../models/divided/vase_antoninaPY.off", "../models/divided/vase_antoninaMY.off"};
    //std::vector<std::string> filenames = {"../models/divided/trophy.off", /*"../models/divided/vase_holeMX.off"*/};
    std::vector<vec3> cameras = {vec3(0, 0, dist), /*vec3(0, 0, -dist),*/ vec3(dist, 0, 0), vec3(-dist, 0, 0), vec3(0, dist, 0), vec3(0, -dist, 0)};
//    std::vector<std::string> filenames = { "../models/analytic/cylinder_lowpoly.off"};
    OctreeNode* root = Fusion::octree_from_samples(openmesh_to_glm(bb_min) - vec3(0.1), octreeSize * 1.1, height,
                                                   filenames, cameras);
    //root = Octree::SimplifyOctree(root, octreeSize/100.0);

    VertexBuffer vertices;
    IndexBuffer indices;
    GenerateMeshFromOctree(/*sphere_octree.*/root, vertices, indices);

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
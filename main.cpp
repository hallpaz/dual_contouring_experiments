#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "Utils.h"

#include "Reconstruction.h"
#include "Contouring.h"
#include "Octree.h"

#include "json.hpp"

using namespace std;
using glm::vec3;

using json = nlohmann::json;
// ----------------------------------------------------------------------------

Real compute_boundingbox(DefaultMesh &mesh, DefaultMesh::Point &bb_min);

const string MESH_FILES = "mesh_files";
const string CAM_POS = "cam_pos";

string threshold_str(Real threshold){
    std::stringstream ss;
    ss << "simp";
    while(threshold < 0){
        ss<< "0";
        threshold *= 10;
    }
    ss << std::to_string(threshold) ;
    return ss.str();
}

int main(int argc, char** argv)
{
    const int height = 6;
//    int dist = 16;

    string folder_name = "../";
    string inputfilename, outputfilename;
    std::cout <<"Input File Name" << endl;
    //std::cin >> inputfilename;
    inputfilename = "ssmbranca.json";
    std::cout <<"Output File Name" << endl;
    std::cin >> outputfilename;

    std::ifstream descriptionStream(folder_name + inputfilename);
    json object_json;
    if (descriptionStream.is_open()) {
        descriptionStream >> object_json;
    }

    std::vector<std::string> filenames = object_json[MESH_FILES];
    std::vector<glm::vec3> cameras;
    for (json::iterator j_it = object_json[CAM_POS].begin(); j_it != object_json[CAM_POS].end(); ++j_it) {
        json cam_json = j_it.value();
        std::cout << cam_json << std::endl;
        cameras.push_back(glm::vec3(cam_json["x"], cam_json["y"], cam_json["z"]));
    }

    DefaultMesh myMesh;
    //OpenMesh::IO::read_mesh(myMesh, "../models/analytic/sphere_lowpoly.off");
    //OpenMesh::IO::read_mesh(myMesh, "../models/divided/vase_antonina.off");
    OpenMesh::IO::read_mesh(myMesh, "../models/branca/k3branca1.off");
    //OpenMesh::IO::read_mesh(myMesh, "../models/taoju/mechanic.off");
    //OpenMesh::IO::read_mesh(myMesh, "../models/cow.off");
    // compute bounding box
    DefaultMesh::Point bb_min;
    Real octreeSize = compute_boundingbox(myMesh, bb_min);


    OctreeNode* root = Fusion::octree_from_samples(openmesh_to_glm(bb_min) - vec3(0.1), octreeSize * 1.1, height,
                                                   filenames, cameras);
    Octree::classify_leaves_vertices(root);

    VertexBuffer vertices;
    IndexBuffer indices;
    GenerateMeshFromOctree(/*sphere_octree.*/root, vertices, indices);


#ifdef DEBUG
    std::cout << "Unoptimized points: " << Octree::unoptimized_points << std::endl;
    std::cout << "Irregular cells " << Octree::irregular_cells << std::endl;
#endif
    std::stringstream filepath;
    filepath << folder_name << outputfilename << "full" << height << ".ply";
    write_Ply(filepath.str(), vertices, indices);


    //generates simplified version
    Real simp_threshold = 100;
    root = Octree::SimplifyOctree(root, simp_threshold);
    GenerateMeshFromOctree(root, vertices, indices);
    std::stringstream simpfilepath;
    simpfilepath << folder_name << threshold_str(simp_threshold) << outputfilename << "full" << height << ".ply";
    write_Ply(simpfilepath.str(), vertices, indices);

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
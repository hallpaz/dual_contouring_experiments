//
// Created by Hallison da Paz on 03/10/2016.
//

#include "Reconstruction.h"
#include "Utils.h"
#include "Constants.h"
#include "old/NormalsEstimator.h"

using glm::vec3;

using std::string;

namespace Fusion
{
    void clean_nodes(OctreeNode* node);

    void clean_nodes(OctreeNode* node)
    {
        if (node == nullptr) return;
        node->clean();
        if (node->type == NODE_INTERNAL){
            for (int i = 0; i < NUM_CHILDREN; ++i) {
                clean_nodes(node->children[i]);
            }
        }
    }


    OctreeNode *octree_from_samples(const glm::vec3 &min, const float size, const unsigned int max_depth,
                                        std::vector<std::string> meshfiles, std::vector<glm::vec3> cameras)
    {
        DefaultMesh mesh;
        OpenMesh::IO::read_mesh(mesh, meshfiles[0]);
        OpenMesh::IO::Options options;
        options += OpenMesh::IO::Options::VertexNormal;
        int i = 0;
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();
        NormalsEstimator::compute_better_normals(mesh);
        /*mesh.request_face_normals();
        mesh.update_face_normals();
        mesh.request_vertex_normals();
        mesh.update_vertex_normals();*/
        OpenMesh::IO::write_mesh(mesh, "../malha_vaso.ply", options);

        Octree demi_octree(min, size, max_depth, mesh, cameras[i++]);
        std::cout << "The first is OK" << std::endl;
        for (std::vector<string>::iterator s_it = meshfiles.begin() + 1; s_it != meshfiles.end(); ++s_it)
        {
            //testing
            //Octree::leafvertexpool.clear();
            // disabled because now we use only intersections

            clean_nodes(demi_octree.root);//TODO: now I can clean the nodes during leaf construction

            DefaultMesh mesh;
            OpenMesh::IO::read_mesh(mesh, *s_it);
            mesh.request_vertex_status();
            mesh.request_edge_status();
            mesh.request_face_status();
            NormalsEstimator::compute_better_normals(mesh);
            /*std::cout << "Opening " << *s_it << std::endl;*/
            Octree::UpdateMeshHierarchy(demi_octree.root, max_depth, mesh);

            //demi_octree.classify_leaves_vertices(cameras[i++], demi_octree.root, mesh);
#ifdef DEBUG
            std::cout << "Divergence: " << Octree::divergence << std::endl;
            std::cout << "Ambiguities solved: " << Octree::ambiguous_vertices << std::endl;
#endif
        }

        return demi_octree.root;
    }

}
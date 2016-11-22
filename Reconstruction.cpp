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

    OctreeNode *update_leaf(OctreeNode *leaf, const DefaultMesh &mesh);
    bool update_children(OctreeNode* node, unsigned int max_depth, const DefaultMesh &mesh);
    OctreeNode* UpdateMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh);
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


    bool update_children(OctreeNode* node, unsigned int max_depth, const DefaultMesh &mesh)
    {
        const float childSize = node->size / 2;
        const int childDepth = node->depth + 1;
        bool hasChildren = false;
        for (int i = 0; i < 8; ++i)
        {
            if (node->children[i] == nullptr)
            {
                OctreeNode* child = new OctreeNode(NODE_INTERNAL,
                                                   node->min + (CHILD_MIN_OFFSETS[i] * childSize),
                                                   childSize, childDepth, node);
                node->children[i] = Octree::BuildMeshHierarchy(child, max_depth, mesh);
            }
            else
            {
                node->children[i] = UpdateMeshHierarchy(node->children[i], max_depth, mesh);
            }
            hasChildren |= (node->children[i] != nullptr);
        }
        return hasChildren;
    }


    OctreeNode* UpdateMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh)
    {
        if (!node) return nullptr;
        //node->clean();
        select_inner_crossing_faces(node, mesh);
        if (node->innerEmpty() && node->crossingEmpty())
        {
            return node;
        }

        if (node->type == NODE_LEAF)
        {
            //won't descend for now, let's see the result in a coarse level
            /*if (node->depth < max_depth && !node->parent->innerEmpty())
            {
                node->type = NODE_INTERNAL;
            }
            else {
                return node;
            }*/
            return node;
        }

        update_children(node, max_depth, mesh);
        return node;
    }


    OctreeNode *octree_from_samples(const glm::vec3 &min, const float size, const unsigned int max_depth,
                                        std::vector<std::string> meshfiles, std::vector<glm::vec3> cameras)
    {
        DefaultMesh mesh;
        OpenMesh::IO::read_mesh(mesh, meshfiles[0]);
        int i = 0;
        /*mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();
        NormalsEstimator::compute_better_normals(mesh);*/

        Octree demi_octree(min, size, max_depth, mesh, cameras[i++]);

        for (std::vector<string>::iterator s_it = meshfiles.begin() + 1; s_it != meshfiles.end(); ++s_it)
        {
            Octree::leafvertexpool.clear();
            clean_nodes(demi_octree.root);
            DefaultMesh mesh;
            OpenMesh::IO::read_mesh(mesh, *s_it);
            /*mesh.request_vertex_status();
            mesh.request_edge_status();
            mesh.request_face_status();
            NormalsEstimator::compute_better_normals(mesh);
            std::cout << "Opening " << *s_it << std::endl;*/
            UpdateMeshHierarchy(demi_octree.root, max_depth, mesh);
            demi_octree.classify_leaves_vertices(cameras[i++], demi_octree.root, mesh);
        }

        return demi_octree.root;
    }

}
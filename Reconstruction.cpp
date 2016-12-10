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

    OctreeNode *update_leaf(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh);
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
        for (int i = 0; i < NUM_CHILDREN; ++i)
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
            if (node->depth < max_depth && !node->parent->innerEmpty())
            {
                node->type = NODE_INTERNAL;
                node->construct_children(max_depth, mesh);
            }
            else {
                node = update_leaf(node, max_depth, mesh);
            }
            return node;
        }

        update_children(node, max_depth, mesh);
        return node;
    }

    OctreeNode *update_leaf(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh)
    {
        // otherwise the voxel contains the surface, so find the edge intersections
        vec3 averageNormal(0.f);
        svd::QefSolver qef;
        bool hasIntersection = false;
        //std::cout << "So far, so good" << std::endl;
        //TODO: optimize computations to avoid redundant intersections (same edge from other cell)
        for (int i = 0; i < NUM_EDGES; ++i) //for each edge
        {
            const vec3 p1 = leaf->get_vertex(edgevmap[i][0]);
            const vec3 p2 = leaf->get_vertex(edgevmap[i][1]);

            vec3 intersection;
            std::vector<vec3> intersection_points, normals;//, face_normals;
            for (std::list<DefaultMesh::FaceHandle>::iterator face = leaf->crossingFaces.begin(); face != leaf->crossingFaces.end(); ++face)
            {
                auto fv_it = mesh.cfv_iter(*face);
                DefaultMesh::VertexHandle a = *fv_it;
                DefaultMesh::VertexHandle b = *(++fv_it);
                DefaultMesh::VertexHandle c = *(++fv_it);

                Vertex face_vertices[3] = {openmesh_to_glm(mesh.point(a)), openmesh_to_glm(mesh.point(b)), openmesh_to_glm(mesh.point(c))};
                //Vertex vertices[3] = { face_vertices[0], face_vertices[1], face_vertices[2]};
                if (moller_triangle_intersection(p1, p2, face_vertices, intersection)) {
                    //keeps the intersection here
                    if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                        continue;
                    }
                    intersection_points.push_back(intersection);

                    float u, v, w;
                    barycentric(intersection, face_vertices[0].position, face_vertices[1].position, face_vertices[2].position, u, v, w);
                    vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                    normals.push_back(glm::normalize(normal_at_intersection));
                }
            }
            if (intersection_points.size() > 1) {
                if (leaf->depth < max_depth){
                    std::cout << intersection_points.size() << " Child Depth: " << leaf->depth+1 << " Child Size: " << leaf->size/2 << std:: endl;

                    if(leaf->construct_children(max_depth, mesh))
                    {
                        leaf->type = NODE_INTERNAL;
                        return leaf;
                    }
                    std::cout << "SERIAO????" << std::endl; //if it has an intersection why not the children?
                    return leaf;
                }
            }
            // if we consider that an intersection happened, we'll consider only the first intersection for now
            if (intersection_points.size() > 0)
            {
                vec3 &n = normals[0];
                vec3 &v = intersection_points[0];
                qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
                averageNormal += n;
                hasIntersection = true;
            }
        }

        //TODO: pass drawinfo data from parent to child when descending
        leaf->drawInfo->qef = leaf->drawInfo->qef + qef.getData();
        leaf->drawInfo->averageNormal += averageNormal;

        leaf->type = NODE_LEAF;
        for (int i = 0; i < NUM_CHILDREN; ++i) {
            std::string vertex_hash = hashvertex(leaf->get_vertex(i));
            if (Octree::leafvertexpool.count(vertex_hash) == 0)
                Octree::leafvertexpool[vertex_hash] = MATERIAL_UNKNOWN;
        }
        return leaf;
    }


    OctreeNode *octree_from_samples(const glm::vec3 &min, const float size, const unsigned int max_depth,
                                        std::vector<std::string> meshfiles, std::vector<glm::vec3> cameras)
    {
        DefaultMesh mesh;
        OpenMesh::IO::read_mesh(mesh, meshfiles[0]);
        int i = 0;
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();
        NormalsEstimator::compute_better_normals(mesh);

        Octree demi_octree(min, size, max_depth, mesh, cameras[i++]);

        std::cout << "The first is OK" << std::endl;
        for (std::vector<string>::iterator s_it = meshfiles.begin() + 1; s_it != meshfiles.end(); ++s_it)
        {
            Octree::leafvertexpool.clear();
            clean_nodes(demi_octree.root);
            DefaultMesh mesh;
            OpenMesh::IO::read_mesh(mesh, *s_it);
            mesh.request_vertex_status();
            mesh.request_edge_status();
            mesh.request_face_status();
            NormalsEstimator::compute_better_normals(mesh);
            /*std::cout << "Opening " << *s_it << std::endl;*/
            UpdateMeshHierarchy(demi_octree.root, max_depth, mesh);
            demi_octree.classify_leaves_vertices(cameras[i++], demi_octree.root, mesh);
        }

        return demi_octree.root;
    }

}
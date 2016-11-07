//
// Created by Hallison da Paz on 03/10/2016.
//

#include <fstream>
#include "Reconstruction.h"
#include "Utils.h"
#include "Constants.h"
#include "NormalsEstimator.h"

using glm::vec3;

using std::string;

namespace Fusion
{

    OctreeNode *update_leaf(OctreeNode *leaf, const DefaultMesh &mesh);
    inline bool update_children(OctreeNode* node, const DefaultMesh &mesh);


    OctreeNode* update_octree(OctreeNode *node, const DefaultMesh &mesh)
    {
        //std::cout << "update_octree" << std::endl;
        if (!node)
        {
            std::cout << "Trying to construct empty node" << std::endl;
            return nullptr;
        }

        select_inner_crossing_faces(node, mesh);
        bool inner_is_empty = node->innerFaces.empty();
        bool crossing_is_empty = node->crossingFaces.empty();
        // cases 2 and 5
        if ( inner_is_empty && crossing_is_empty)
        {   //no data here to update. other frame must've created this node
            return clean_node(node);
        }

        bool hasChildren = false;
        if (node->type == NODE_LEAF)
        {

            if (node->parent->innerFaces.empty() || node->height == 0)
            {
                return update_leaf(node, mesh);
            }
            else
            {
                trace("leaf else");
                //return  update_leaf(node, mesh);
                hasChildren = construct_children(node, mesh);
                if (hasChildren)
                {
                    node->type = NODE_INTERNAL;
                }
                return clean_node(node);
                //return update_octree(node, mesh);
            }



            //return update_leaf(node, mesh);



//            if (node->height > 0 && ! node->parent->innerFaces.empty())
//            {
//                hasChildren = construct_children(node, mesh);
//                if (hasChildren)
//                {
//                    node->type = NODE_INTERNAL;
//                }
//                return clean_node(node);
//            }
//            else
//            {
//                if (node->parent->innerFaces.empty() || !crossing_is_empty)
//                {
//                    return update_leaf(node, mesh);
//                }
//                return clean_node(node);
//            }
            // case 3
            if (crossing_is_empty) {
                //if (node->height > 0)
                //{
                    hasChildren = construct_children(node, mesh);
                    if (hasChildren)
                    {
                        node->type = NODE_INTERNAL;
                    }
                //}
                return clean_node(node);
            }
            //case 4
            if (!crossing_is_empty && inner_is_empty) {
                if (/*node->height > 0 && */ !node->parent->innerFaces.empty())
                {
                    hasChildren = construct_children(node, mesh);
                    if (hasChildren)
                    {
                        node->type = NODE_INTERNAL;
                    }
                    //clear memory used for inner and crossing faces
                    return clean_node(node);
                }
                return update_leaf(node, mesh);
            }
            //case 6
            if (!crossing_is_empty && !inner_is_empty) {
                //if (node->height > 0) {
                    hasChildren = construct_children(node, mesh);
                    if (hasChildren)
                    {
                        node->type = NODE_INTERNAL;
                        return clean_node(node);
                    }//ATTENTION
                    //clear memory used for inner and crossing faces

                //}
                trace("CARILHO!!!!!!");
                return update_leaf(node, mesh);
            }
        }
        //case 1
        //trace("caso 1: update_children direto");
        update_children(node, mesh);
        return clean_node(node);
    }

    int compute_test_sign(vec3 cube_vertex)
    {
        int sign = glm::length(cube_vertex) < 8.0f ? MATERIAL_SOLID : MATERIAL_AIR;
        return sign;
    }

    OctreeNode *update_leaf(OctreeNode *leaf, const DefaultMesh &mesh) {
        if (!leaf)
        {
            std::cout << "Trying to construct a leaf in the middle" << std::endl;
            return nullptr;
        }
        //trace("update_leaf");
        // otherwise the voxel contains the surface, so find the edge intersections
        vec3 averageNormal(0.f);
        svd::QefSolver qef;
        bool hasIntersection = false;
        //int corners = leaf->drawInfo->corners;
        int corners = 0;

        //vertices classification
        //TODO: first we'll try to "merge" both signs using logical operations. if it doens't work, we'll have to figure out a way to update sign by sign
        int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                           MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
        for (int i = 0; i < 12; ++i) //for each edge
        {
            const int c1 = edgevmap[i][0];
            const int c2 = edgevmap[i][1];
            const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
            const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);
            int sign1 = compute_test_sign(p1);
            int sign2 = compute_test_sign(p2);

            vecsigns[c1] = sign1; vecsigns[c2] = sign2;

            if (sign1 == sign2)
            {
                continue;
            }

            vec3 intersection;
            std::vector<vec3> intersection_points, normals, face_normals;
            for (std::list<DefaultMesh::FaceHandle>::iterator face = leaf->crossingFaces.begin(); face != leaf->crossingFaces.end(); ++face)
            {
                auto fv_it = mesh.cfv_iter(*face);
                DefaultMesh::VertexHandle a = *fv_it;
                DefaultMesh::VertexHandle b = *(++fv_it);
                DefaultMesh::VertexHandle c = *(++fv_it);

                vec3 face_vertices[3] = {openmesh_to_glm(mesh.point(a)), openmesh_to_glm(mesh.point(b)), openmesh_to_glm(mesh.point(c))};
                Vertex vertices[3] = { face_vertices[0], face_vertices[1], face_vertices[2]};

                if (moller_triangle_intersection(p1, p2, vertices, intersection)) {
                    //keeps the intersection here
                    if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                        continue;
                    }
                    intersection_points.push_back(intersection);

                    float u, v, w;
                    barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);

                    vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                    normal_at_intersection =  glm::normalize(normal_at_intersection);
                    normals.push_back(normal_at_intersection);
                    vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                    face_normals.push_back(face_normal);
                }
            }

            for (int j = 0; j < 8; ++j) {
                const vec3 cube_vertex = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[j]);
                vecsigns[j] = compute_test_sign(cube_vertex);//glm::length(cube_vertex) < 8.0f ? MATERIAL_SOLID : MATERIAL_AIR;
            }

            // if we consider that an intersection happened.
            if ((intersection_points.size() > 0) && (vecsigns[c1] != vecsigns[c2]))
            {   // we'll consider only the first intersection for now
                vec3 &n = normals[0];
                vec3 &v = intersection_points[0];
                qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
                averageNormal += n;
                hasIntersection = true;
            }
        }

//        updateSignsArray(vecsigns, 8);

        for (size_t i = 0; i < 8; ++i)
        {   //encode the signs to the corners variable to save memory
            corners |= (vecsigns[i] << i);
  //          updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
        }

        /*if (corners == 0 || corners == 255)
            hasIntersection = false;
        else
            hasIntersection = true;*/

        if (!hasIntersection)
        {   // can't delete leaf here, because other image constructed it
            if (corners != 0 && corners != 255){
                trace("Update should have intersected");
                auto p = leaf->min + vec3(leaf->size/2, leaf->size/2, leaf->size/2);
                auto n = vec3(0, 0, 1);
                qef.add(p.x, p.y, p.z, n.x, n.y, n.z);
                averageNormal = n;
                hasIntersection = true;
            }
            else {
                return clean_node(leaf);
            }
        }

        if (leaf->drawInfo->corners != corners){
            std::ofstream signfile("../signs_bug.txt");
            if (signfile.is_open()){
                for (int i = 0; i < 8; ++i) {
                    signfile << ((leaf->drawInfo->corners >> i) & 1);
                }
                signfile << " ";
                for (int i = 0; i < 8; ++i) {
                    signfile << ((corners >> i) & 1);
                }
                signfile << std::endl;
            }
            else {
                std::cout << "XABU SIGNFILE" << std::endl;
            }
            signfile.close();
        }

        leaf->drawInfo->averageNormal += averageNormal;


        leaf->drawInfo->corners |= corners;
        //leaf->drawInfo->corners = corners;

        leaf->drawInfo->qef.add(qef.getData());
        //leaf->drawInfo->qef = qef.getData();

        return clean_node(leaf);
    }

    OctreeNode* octree_from_samples(const glm::vec3 &min, const float size, const int height, std::vector<string> meshfiles)
    {
        DefaultMesh mesh;
        OpenMesh::IO::read_mesh(mesh, meshfiles[0]);
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();
        NormalsEstimator::compute_better_normals(mesh);

        OctreeNode* root = BuildOctreeFromOpenMesh(min, size, height, mesh);
        float simpthreshold = size/100;
        //root = SimplifyOctree(root, simpthreshold);
        for (std::vector<string>::iterator s_it = meshfiles.begin() + 1; s_it != meshfiles.end(); ++s_it)
        {
            DefaultMesh mesh;
            OpenMesh::IO::read_mesh(mesh, *s_it);
            mesh.request_vertex_status();
            mesh.request_edge_status();
            mesh.request_face_status();
            NormalsEstimator::compute_better_normals(mesh);
            std::cout << "Opening " << *s_it << std::endl;
            root = update_octree(root, mesh);
            // we must clear the data used to build the representation from the previous mesh
            OctreeNode::vertexpool.clear();
            OctreeNode::edgepool.clear();
            //root = SimplifyOctree(root, simpthreshold);
        }

        return root;
    }

    inline bool update_children(OctreeNode* node, const DefaultMesh &mesh)
    {
        const float childSize = node->size / 2;
        const int childHeight = node->height - 1;
        bool hasChildren = false;
        for (int i = 0; i < 8; ++i)
        {
            if (node->children[i] == nullptr)
            {
                OctreeNode* child = new OctreeNode(NODE_INTERNAL,
                                                   node->min + (CHILD_MIN_OFFSETS[i] * childSize),
                                                   childSize, childHeight, node);
                node->children[i] = ConstructOctreeNodesFromOpenMesh(child, mesh);
            }
            else
            {
                node->children[i] = update_octree(node->children[i], mesh);
            }
            hasChildren |= (node->children[i] != nullptr);
        }
        return hasChildren;
    }

}


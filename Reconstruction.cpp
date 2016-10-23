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
                if (node->height > 0)
                {
                    hasChildren = construct_children(node, mesh);
                    if (hasChildren)
                    {
                        node->type = NODE_INTERNAL;
                    }
                }
                return clean_node(node);
            }
            //case 4
            if (!crossing_is_empty && inner_is_empty) {
                if (node->height > 0 && ! node->parent->innerFaces.empty())
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
                if (node->height > 0) {
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
        }
        //case 1
        update_children(node, mesh);
        return clean_node(node);
    }

    OctreeNode *update_leaf(OctreeNode *leaf, const DefaultMesh &mesh) {
        if (!leaf)
        {
            std::cout << "Trying to construct a leaf in the middle" << std::endl;
            return nullptr;
        }
        trace("update_leaf");
        //std::cout << "Leaf height: " << leaf->height << std::endl;
        // otherwise the voxel contains the surface, so find the edge intersections
        vec3 averageNormal(0.f);
        svd::QefSolver qef;
        bool hasIntersection = false;
        //int corners = leaf->drawInfo->corners;
        int corners = 0;
        //int vecsigns[8];
        /*for (int j = 0; j < 8; ++j) { //we'll begin from the already known corner
            vecsigns[j] = (corners >> j) & 1;
        }*/
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

            //computes hash for edge
            std::string edgehash = hashedge(p1, p2);
            HermiteData edgedata;
            if (OctreeNode::edgepool.count(edgehash) != 0)
            {
                // if edge data already exists, retrieve the data
                edgedata = OctreeNode::edgepool[edgehash];
                if (edgedata.hasIntersection()){
                    vecsigns[c1] = OctreeNode::vertexpool[hashvertex(p1)];
                    vecsigns[c2] = OctreeNode::vertexpool[hashvertex(p2)];
                    if (vecsigns[c1] == MATERIAL_UNKNOWN || vecsigns[c2] == MATERIAL_UNKNOWN)
                    {
                        std::cout << "SIGNAL INCONSISTENCY FOR VERTICES IN INTERSECTION" << std::endl;
                    }
                }
                //continue;
            }

            vec3 intersection;
            std::vector<vec3> intersection_points;
            std::vector<vec3> normals;
            std::vector<vec3> face_normals;
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

                    if (vecsigns[c1] != MATERIAL_SOLID){
                        vecsigns[c1] = computeSideOfPoint(p1, intersection, face_normal);
                    }
                    if (vecsigns[c2] != MATERIAL_SOLID){
                        vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
                    }
                    /*vecsigns[c1] = glm::length(p1) < 8.0f ? MATERIAL_SOLID : MATERIAL_AIR;
                    vecsigns[c2] = glm::length(p2) < 8.0f ? MATERIAL_SOLID : MATERIAL_AIR;*/
                }
            }
            if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
                if (intersection_points.size()%2 == 0)
                {
                    const float childSize = leaf->size / 2;
                    const int childHeight = leaf->height - 1;
                    leaf->type = NODE_INTERNAL;
                    std::cout << intersection_points.size() << " Child Height: " << childHeight << " Child Size: " << childSize << std:: endl;
                    bool hasChildren = construct_children(leaf, mesh);
                    if (!hasChildren)
                    {
                        delete leaf;
                        return nullptr;
                    }

                    leaf->crossingFaces.clear();
                    leaf->innerFaces.clear();
                    return leaf;
                }
                else
                {
                    // they are on the opposite side of the surface
                    int nearindex = glm::distance(p1, intersection_points[0]) < glm::distance(p1, intersection_points[1]) ? 0 : 1;
                    nearindex = glm::distance(p1, intersection_points[nearindex]) < glm::distance(p1, intersection_points[2]) ? nearindex : 2;

                    //vecsigns[c1] = computeSideOfPoint(p1, intersection_points[nearindex], face_normals[nearindex]);
                    // they are on the same side of the surface. We must ignore the intersection
                    //vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
                    std::cout << "CHEGA MAIS" << std::endl;
                }
            }


            for (int j = 0; j < 8; ++j) {
                const vec3 cube_vertex = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[j]);
                vecsigns[j] = glm::length(cube_vertex) < 8.0f ? MATERIAL_SOLID : MATERIAL_AIR;
            }


            // if we consider that an intersection happened.
            if ((intersection_points.size() > 0) && (vecsigns[c1] != vecsigns[c2]))
            {   // we'll consider only the first intersection for now
                vec3 &n = normals[0];
                vec3 &v = intersection_points[0];
                qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
                averageNormal += n;
                hasIntersection = true;
                edgedata.intersection = intersection_points[0];
                edgedata.normal = normals[0];
            }

            OctreeNode::edgepool[edgehash] = edgedata;
        }

        if (!hasIntersection)
        {   // can't delete leaf here, because other image constructed it
            return clean_node(leaf);
        }
//        updateSignsArray(vecsigns, 8);

        for (size_t i = 0; i < 8; ++i)
        {   //encode the signs to the corners variable to save memory
            corners |= (vecsigns[i] << i);
  //          updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
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
        //root = SimplifyOctree(root, size/10);
        float simpthreshold = size/100;
        for (std::vector<string>::iterator s_it = meshfiles.begin() + 1; s_it != meshfiles.end(); ++s_it)
        {
            //root = SimplifyOctree(root, simpthreshold);
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
            //root = SimplifyOctree(root, size/10);
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


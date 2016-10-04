//
// Created by Hallison da Paz on 03/10/2016.
//

#include "Reconstruction.h"
#include "Utils.h"
#include "Constants.h"
#include "NormalsEstimator.h"

using glm::vec3;

using std::string;

namespace Fusion
{

    OctreeNode *update_leaf(OctreeNode *leaf, const DefaultMesh &mesh);


    OctreeNode* update_octree(OctreeNode *node, const DefaultMesh &mesh)
    {
        std::cout << "update_octree" << std::endl;
        if (!node)
        {
            std::cout << "Trying to construct empty node" << std::endl;
            return nullptr;
        }

        select_inner_crossing_faces(node, mesh);

        if (node->innerFaces.size() == 0)
        {   //Empty space, no triangles crossing or inside this cell
            if (node->crossingFaces.size() == 0)
            {
                //we should only delete this node if it didn't exist before. If it existed, either it has children or it is a leaf
                bool hasChildren = false;
                for (int i = 0; i < 8; ++i) {
                    if (node->children[i] != nullptr)
                    {
                        hasChildren = true;
                    }
                }
                if (hasChildren || node->type == Node_Leaf)
                    return node;

                delete node;
                return nullptr;
            }
            else
            {
                if (node->parent->innerFaces.size() == 0 && node->type == Node_Leaf)
                {
                    std::cout << "going to update_leaf" << std::endl;
                    return update_leaf(node, mesh);
                }

                //if it wasn't a leaf, we let it unchanged
                std::cout << "no changes" << std::endl;
                return node;
            }
        }

        if (node->height == 0)
        {
            //std::cout << "HEIGHT 0 ACTIVATED!!" << std::endl;
            //this is a moment to construct, not update
            if (node->type == Node_Leaf){
                std::cout << "Heigh 0, will update" << std::endl;
                return update_leaf(node, mesh);
            }

            std::cout << "Height 0, will construct" << std::endl;
            return ConstructLeafFromOpenMesh(node, mesh);
        }

        const float childSize = node->size / 2;
        const int childHeight = node->height - 1;
        bool hasChildren = false;

        for (int i = 0; i < 8; i++)
        {
            if (node->children[i] == nullptr)
            {
                OctreeNode* child = new OctreeNode;
                child->size = childSize;
                child->height = childHeight;
                child->min = node->min + (CHILD_MIN_OFFSETS[i] * childSize);
                child->type = Node_Internal;
                child->parent = node;

                std::cout << i << " constructing inside update " << node->type << "Inner " << node->innerFaces.size() << "Crossing: " << node->crossingFaces.size() << std::endl;
                node->children[i] = ConstructOctreeNodesFromOpenMesh(child, mesh);
            }
            else
            {
                update_octree(node->children[i], mesh);
            }
            hasChildren |= (node->children[i] != nullptr);
        }

        if (!hasChildren)
        {
            std::cout << "will delete node" << std::endl;
            delete node;
            return nullptr;
        }

        node->crossingFaces.clear();
        node->innerFaces.clear();

        std::cout << "Node updated" << std::endl;
        return node;
    }

    OctreeNode *update_leaf(OctreeNode *leaf, const DefaultMesh &mesh) {
        if (!leaf)
        {
            std::cout << "Trying to construct a leaf in the middle" << std::endl;
            return nullptr;
        }
        //std::cout << "Leaf height: " << leaf->height << std::endl;
        // otherwise the voxel contains the surface, so find the edge intersections
        vec3 averageNormal(0.f);
        svd::QefSolver qef;
        bool hasIntersection = false;
        int corners = 0;
        //vertices classification
        //TODO: first we'll try to "merge" both signs using logical operations. if it doens't work, we'll have to figure out a way to update sign by sign
        int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                           MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};

        for (int i = 0; i < 12; i++) //for each edge
        {
            const int c1 = edgevmap[i][0];
            const int c2 = edgevmap[i][1];
            const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
            const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

            //computes hash for edge
            std::string edgehash = hashedge(p1, p2);
            HermiteData edgedata;
            if (OctreeNode::edgepool.count(edgehash) != 0){
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
            std::cout << "Update Leaf" << std::endl;
            for (std::list<DefaultMesh::FaceHandle>::iterator face = leaf->crossingFaces.begin(); face != leaf->crossingFaces.end(); ++face)
            {
                if (! mesh.is_valid_handle(*face))
                {
                    std::cout << "Invalid Handle" << std::endl;
                    continue;
                }
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

                    if (vecsigns[c1] != MATERIAL_SOLID/*vecsigns[c1] == MATERIAL_UNKNOWN*/){
                        vecsigns[c1] = computeSideOfPoint(p1, intersection, face_normal);
                    }
                    if (vecsigns[c2] != MATERIAL_SOLID/*vecsigns[c2] == MATERIAL_UNKNOWN*/){
                        vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
                    }
                }
            }
            if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
                if (intersection_points.size()%2 == 0){
                    const float childSize = leaf->size / 2;
                    const int childHeight = leaf->height - 1;
                    std::cout << intersection_points.size() << " Child Height: " << childHeight << " Child Size: " << childSize << std:: endl;
                    bool hasChildren = false;

                    for (int i = 0; i < 8; i++)
                    {
                        OctreeNode* child = new OctreeNode;
                        child->size = childSize;
                        child->height = childHeight;
                        child->min = leaf->min + (CHILD_MIN_OFFSETS[i] * childSize);
                        child->type = Node_Internal;
                        child->parent = leaf;

                        leaf->children[i] = ConstructOctreeNodesFromOpenMesh(child, mesh);
                        hasChildren |= (leaf->children[i] != nullptr);
                    }

                    if (!hasChildren)
                    {
                        delete leaf;
                        return nullptr;
                    }

                    leaf->crossingFaces.clear();
                    leaf->innerFaces.clear();
                    return leaf;
                }
                else{
                    // they are on the opposite side of the surface
                    int nearindex = glm::distance(p1, intersection_points[0]) < glm::distance(p1, intersection_points[1]) ? 0 : 1;
                    nearindex = glm::distance(p1, intersection_points[nearindex]) < glm::distance(p1, intersection_points[2]) ? nearindex : 2;
                    vecsigns[c1] = computeSideOfPoint(p1, intersection_points[nearindex], face_normals[nearindex]);
                    // they are on the same side of the surface. We must ignore the intersection
                    vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;

                    std::cout << "CHEGA MAIS" << std::endl;
                }


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
        {   // can't delet leaf here, because other image constructed it
            leaf->crossingFaces.clear();
            leaf->innerFaces.clear();
            return leaf;
        }
        updateSignsArray(vecsigns, 8);

        for (size_t i = 0; i < 8; i++)
        {   //encode the signs to the corners variable to save memory
            corners |= (vecsigns[i] << i);
            updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
        }


        leaf->drawInfo->averageNormal += averageNormal;
        leaf->drawInfo->corners |= corners;

        // redundant??
        leaf->type = Node_Leaf;

        leaf->crossingFaces.clear();
        leaf->innerFaces.clear();

        return leaf;
    }

    OctreeNode* octree_from_samples(const glm::vec3 &min, const float size, const int height, std::vector<string> meshfiles)
    {

        OctreeNode* root = new OctreeNode();
        root->min = min;
        root->size = size;
        root->height = height;
        root->type = Node_Internal;

        DefaultMesh mesh;
        OpenMesh::IO::read_mesh(mesh, meshfiles[0]);
//        OpenMesh::IO::read_mesh(mesh, "../models/sphere8.off");
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();
        NormalsEstimator::compute_better_normals(mesh);

        root = BuildOctreeFromOpenMesh(min, size, height, mesh);
        root = SimplifyOctree(root, size/10);
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
            root = SimplifyOctree(root, size/10);
        }


        return root;
    }

}


//
// Created by Hallison da Paz on 03/10/2016.
//

#include "Reconstruction.h"
#include "Utils.h"
#include "Constants.h"

using glm::vec3;

using std::string;

namespace Fusion
{

    OctreeNode *ConstructLeafFusion(OctreeNode *leaf, const DefaultMesh &mesh);


    void divideFacesByLocation(OctreeNode *node, std::list<DefaultMesh::FaceHandle> &facesList, const DefaultMesh &mesh)
    {
        auto f_it = facesList.begin();
        //parent's inner triangles
        while(f_it != facesList.end())
        {
            switch (triangleRelativePosition(mesh, *f_it, node->min, node->size))
            {
                case INSIDE:
                    node->innerFaces.push_back(*f_it);
                    /*If the triangle is located inside the cell, we remove it from the cell's parent list*/
                    f_it = facesList.erase(f_it);
                    break;
                case CROSSING:
                    node->crossingFaces.push_back(*f_it);
                    /*If the triangle might cross the cell, we can't remove it from the cell's parent list
                     * because it might cross other cells as well*/
                    ++f_it;
                    break;
                default:
                    ++f_it;
            }
        }
    }


    OctreeNode* update_octree(OctreeNode *node, const DefaultMesh &mesh)
    {
        if (!node)
        {
            std::cout << "Trying to construct empty node" << std::endl;
            return nullptr;
        }

        if (node->parent != nullptr)
        {
            divideFacesByLocation(node, node->parent->innerFaces, mesh); //parent's inner triangles
            divideFacesByLocation(node, node->parent->crossingFaces, mesh); //parent's crossing triangles
        }
        else
        {   //initializes the parent list with all triangles
            for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
            {
                node->innerFaces.push_back(*f_it);
            }
        }

        if (node->innerFaces.size() == 0)
        {   //Empty space, no triangles crossing or inside this cell
            if (node->crossingFaces.size() == 0)
            {
                //std::cout << "Empty space HERE!!" << std::endl;
                delete node;
                return nullptr;
            }
            else
            {
                if (node->parent->innerFaces.size() == 0)
                    return ConstructLeafFusion(node, mesh);
            }
        }

        if (node->height == 0)
        {
            //std::cout << "HEIGHT 0 ACTIVATED!!" << std::endl;
            return ConstructLeafFusion(node, mesh);
        }

        const float childSize = node->size / 2;
        const int childHeight = node->height - 1;
        bool hasChildren = false;

        for (int i = 0; i < 8; i++)
        {
            OctreeNode* child = new OctreeNode;
            child->size = childSize;
            child->height = childHeight;
            child->min = node->min + (CHILD_MIN_OFFSETS[i] * childSize);
            child->type = Node_Internal;
            child->parent = node;

            node->children[i] = update_octree(child, mesh);
            hasChildren |= (node->children[i] != nullptr);
        }

        if (!hasChildren)
        {
            delete node;
            return nullptr;
        }

        return node;
    }

    OctreeNode *ConstructLeafFusion(OctreeNode *leaf, const DefaultMesh &mesh) {
        if (!leaf)
        {
            std::cout << "Trying to construct a leaf in the middle" << std::endl;
            return nullptr;
        }
        //std::cout << "Leaf height: " << leaf->height << std::endl;
        // otherwise the voxel contains the surface, so find the edge intersections
        vec3 averageNormal(0.f);
        svd::QefSolver qef;
        svd::QefSolver featureQef;
        bool hasIntersection = false;
        int corners = 0;
        //vertices classification
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
//                averageNormal += edgedata.normal;
//                qef.add(edgedata.intersection.x, edgedata.intersection.y, edgedata.intersection.z,
//                        edgedata.normal.x, edgedata.normal.y, edgedata.normal.z);
//                hasIntersection = true;
                    vecsigns[c1] = OctreeNode::vertexpool[hashvertex(p1)];
                    vecsigns[c2] = OctreeNode::vertexpool[hashvertex(p2)];
                    if (vecsigns[c1] == MATERIAL_UNKNOWN || vecsigns[c2] == MATERIAL_UNKNOWN)
                    {
                        std::cout << "SIGNAL INCONSISTENCY FOR VERTICES IN INTERSECTION" << std::endl;
                    }
//                vecsigns[c1] = computeSideOfPoint(p1, edgedata.intersection, edgedata.normal);
//                vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
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
//                int nearindex = glm::distance(p1, intersection_points[0]) < glm::distance(p1, intersection_points[1]) ? 0 : 1;
//                vecsigns[c1] = computeSideOfPoint(p1, intersection_points[nearindex], face_normals[nearindex]);
//                // they are on the same side of the surface. We must ignore the intersection
//                vecsigns[c2] = vecsigns[c1];
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
        {   // voxel is full inside or outside the volume
            delete leaf;
            return nullptr;
        }
        updateSignsArray(vecsigns, 8);

        for (size_t i = 0; i < 8; i++)
        {   //encode the signs to the corners variable to save memory
            corners |= (vecsigns[i] << i);
            updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
        }

//        std::ofstream interiorfile, exteriorfile;
//        //gridfile.open("/home/hallpaz/Workspace/dual_contouring_experiments/grid_color_updated.ply", std::ios::app);
//        interiorfile.open("../subproducts/interior_color_updated.ply", std::ios::app);
//        exteriorfile.open("../subproducts/exterior_color_updated.ply", std::ios::app);
//        for (size_t i = 0; i < 8; i++) {
//            corners |= (vecsigns[i] << i);
//
//            const vec3 cornerPos = leaf->min + CHILD_MIN_OFFSETS[i]*leaf->size;
//            if (vecsigns[i] == MATERIAL_SOLID) {
//                interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 128 << " " << 255 << std::endl;
//            }
//
//            if (vecsigns[i] == MATERIAL_AIR) {
//                exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 64 << " " << 255 << " " << 64 << std::endl;
//            }
//
//            if (vecsigns[i] == MATERIAL_UNKNOWN) {
//                interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
//                exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
//            }
//        }
//        interiorfile.close();
//        exteriorfile.close();



        //qef.setData(qef.getData()*0.4f + featureQef.getData()*0.6f);
        //qef.add(featureQef.getData());



//        svd::Vec3 qefPosition;
//        qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
//
//        OctreeDrawInfo* drawInfo = new OctreeDrawInfo;
//        drawInfo->position = vec3(qefPosition.x, qefPosition.y, qefPosition.z);
//        drawInfo->qef = qef.getData();
//
//        const vec3 min = vec3(leaf->min);
//        const vec3 max = vec3(leaf->min + vec3(leaf->size));
//        if (drawInfo->position.x < min.x || drawInfo->position.x > max.x ||
//            drawInfo->position.y < min.y || drawInfo->position.y > max.y ||
//            drawInfo->position.z < min.z || drawInfo->position.z > max.z)
//        {
//            const auto& mp = qef.getMassPoint();
//            drawInfo->position = vec3(mp.x, mp.y, mp.z);
//        }




        leaf->drawInfo->averageNormal += averageNormal;
        leaf->drawInfo->corners |= corners;

        leaf->type = Node_Leaf;
        //leaf->drawInfo = drawInfo;

        return leaf;
    }

    OctreeNode* octree_from_samples(const glm::vec3 &min, const float size, const int height, std::vector<string> meshfiles)
    {

        OctreeNode* root = new OctreeNode();
        root->min = min;
        root->size = size;
        root->height = height;
        root->type = Node_Internal;



        for (std::vector<string>::iterator s_it = meshfiles.begin(); s_it != meshfiles.end(); ++s_it)
        {
            DefaultMesh mesh;
            OpenMesh::IO::read_mesh(mesh, *s_it);
            std::cout << "Opening " << *s_it << std::endl;
            update_octree(root, mesh);
        }



        return root;
    }

}


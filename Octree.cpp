//
// Created by Hallison da Paz on 18/11/2016.
//
/*

Implementations of Octree member functions.

Copyright (C) 2011  Tao Ju

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
(LGPL) as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <fstream>
#include "Octree.h"
#include "Utils.h"


using glm::vec3;
// ----------------------------------------------------------------------------
std::unordered_map<std::string, int> Octree::leafvertexpool;
std::unordered_map<std::string, HermiteData> Octree::edgepool;

#ifdef DEBUG
int Octree::unoptimized_points = 0;
int Octree::divergence = 0;
int Octree::ambiguous_vertices = 0;
int Octree::irregular_cells = 0;
#endif

// ----------------------------------------------------------------------------
int classify_vertex(glm::vec3 vertex, glm::vec3 cam_origin, OctreeNode* root, DefaultMesh &mesh);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
bool is_good_optimization(OctreeNode *node, svd::QefData qefdata)
{
    svd::Vec3 qefPosition;
    svd::QefSolver qef;
    qef.setData(qefdata);
    if (qefdata.numPoints > 0)
    {
        qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
        vec3 optmized_position = vec3(qefPosition.x, qefPosition.y, qefPosition.z);
        const vec3 min = vec3(node->min);
        const vec3 max = vec3(node->min + vec3(node->size));
        if (optmized_position.x < min.x || optmized_position.x > max.x ||
            optmized_position.y < min.y || optmized_position.y > max.y ||
            optmized_position.z < min.z || optmized_position.z > max.z)
        {
            return false;
        }
    } else{
        return false;
    }
    return true;
}


bool OctreeNode::construct_or_update_children(unsigned int max_depth, const DefaultMesh &mesh)
{
    const float childSize = this->size / 2;
    const int childHeight = this->depth + 1;
    bool hasChildren = false;
    for (int i = 0; i < NUM_CHILDREN; i++)
    {
        if (this->children[i] == nullptr) {
            OctreeNode *child = new OctreeNode(NODE_INTERNAL,
                                               this->min + (CHILD_MIN_OFFSETS[i] * childSize),
                                               childSize,
                                               childHeight,
                                               this);
// TODO: check if this can help - seems to not
//            if (this->drawInfo)
//            {
//                child->drawInfo = new OctreeDrawInfo();
//                child->drawInfo->qef = child->drawInfo->qef + this->drawInfo->qef;
//                child->drawInfo->averageNormal += this->drawInfo->averageNormal;
//            }
            this->children[i] = Octree::BuildMeshHierarchy(child, max_depth, mesh);
        }
        else {
            this->children[i] = Octree::UpdateMeshHierarchy(this->children[i], max_depth, mesh);
        }
        hasChildren |= (this->children[i] != nullptr);
    }
    return hasChildren;
}

Octree::Octree(glm::vec3 min, Real size, unsigned int max_depth, DefaultMesh &mesh, vec3 cam_origin)
{
    root = new OctreeNode(NODE_INTERNAL, min, size, 0);
    trace("building hierarchy");
    BuildMeshHierarchy(root, max_depth, mesh);

#ifdef DEBUG
    std::cout << "Divergence: " << Octree::divergence << std::endl;
    std::cout << "Ambiguities solved: " << Octree::ambiguous_vertices << std::endl;
#endif
}

OctreeNode *Octree::BuildMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh)
{
    if (!node) return nullptr;

    select_inner_crossing_faces(node, mesh);
    if (node->innerEmpty() && node->crossingEmpty())
    {   //Empty space, no triangles crossing or inside this cell
        delete node;
        return nullptr;
    }
    if ((/*(node->parent && node->parent->innerEmpty()) ||*/ node->depth == max_depth) && node->depth > 4)
    {
        return construct_or_update_leaf(node, max_depth, mesh);
    }
    if (node->construct_or_update_children(max_depth, mesh))
    {
        return node;
    }

    delete node;
    return nullptr;
}

OctreeNode* Octree::UpdateMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh)
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
        if (node->depth < max_depth /*&& !node->parent->innerEmpty()*/)
        {
            node->type = NODE_INTERNAL;
            node->construct_or_update_children(max_depth, mesh);
        }
        else
        {
            node = construct_or_update_leaf(node, max_depth, mesh);
        }
        return node;
    }

    node->construct_or_update_children(max_depth, mesh);
    return node;
}

OctreeNode *Octree::construct_or_update_leaf(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh)
{
    if (leaf == nullptr)
        return nullptr;
    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);

    svd::QefSolver qef;
    svd::QefSolver tmpQef;

    bool hasIntersection = false;
    int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                       MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};

    //TODO: optimize computations to avoid redundant intersections (same edge from other cell)
    int edges_intersected = 0;

    std::cout << "LET THE GAME BEGIN!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl << std::endl;
    bool inversion_flag = false;
    for (int i = 0; i < NUM_EDGES; ++i) //for each edge
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

        // ----------------------------------------------------------------------------
        //Check if this edge was already evaluated
        /*std::string edgehash = hashedge(p1, p2);
        if (edgepool.count(edgehash) != 0){
            vecsigns[c1] = leafvertexpool[hashvertex(p1)];
            vecsigns[c2] = leafvertexpool[hashvertex(p2)];
            HermiteData data = edgepool[edgehash];
            if (data.is_valid()){
                auto v = data.intersection;
                auto n = data.normal;
                qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
                averageNormal += n;
                hasIntersection = true;
                continue;
            }

        }*/
        // End check
        // ----------------------------------------------------------------------------

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
            //trace("intersection");
            if (moller_triangle_intersection(p1, p2, vertices, intersection)) {
                //keeps the intersection here
                if (intersection_points.size() > 0)
                {
                    bool redundant_intersection = false;
                    for(int i = 0; i < intersection_points.size(); ++i)
                    {
                        if (glm::distance(intersection, intersection_points[i]) < POINT_DISTANCE_THRESHOLD)
                        {
                            redundant_intersection = true;
                        }
                    }
                    if(redundant_intersection)
                    {
                        continue;
                    }
                }
                intersection_points.push_back(intersection);

                float u, v, w;
                if(u < 0){
                    u = 0.0;
                    leaf->dangerous_intersection = true;
                }
                barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);
                //u = 1.0f; v = 1.0f; w = 1.0f;
                vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                normal_at_intersection =  glm::normalize(normal_at_intersection);
                normals.push_back(normal_at_intersection);
                vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                face_normals.push_back(face_normal);
            }
        }
        if (intersection_points.size() > 1)
        {
            if (leaf->depth < max_depth)
            {
                if(leaf->construct_or_update_children(max_depth, mesh))
                {
                    leaf->type = NODE_INTERNAL;
                    return leaf;
                }
                std::cout << "SERIAO????" << std::endl; //if it has an intersection why not the children?
                exit(14);
                if (leaf->type == NODE_LEAF) //it's being updated
                {
                    return leaf;
                }
                delete leaf;
                return nullptr;
            }
        }
        // if we consider that an intersection happened.
        // we'll consider only the first intersection for now
        if (intersection_points.size() > 0)
        {
            std::cout << "Intersected Edge: " << c1 << "--" << c2 << std::endl;
            bool should_add = true;
            /*BEGIN VERTEX CLASSIFICATION*/
            int sign1 = MATERIAL_UNKNOWN;
            int sign2 = MATERIAL_UNKNOWN;

            if (intersection_points.size() > 1)
            {

                int nearmost_index1 = compute_nearmost_index(p1, intersection_points);
                sign1 = computeSideOfPoint(p1, intersection_points[nearmost_index1], face_normals[nearmost_index1]);

                if (intersection_points.size()%2 == 0){
                    sign2 = sign1;
                }
                else
                {
                    sign2 = sign1 == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
                }

                // PROVA REAL
                int nearmost_index2 = compute_nearmost_index(p2, intersection_points);
                int corrector = computeSideOfPoint(p2, intersection_points[nearmost_index2], face_normals[nearmost_index2]);

                if (corrector != sign2){
                    std::cout << "DIVERGENCIA NO CALCULO DOS SINAIS EM CASO MULTIPLO" << std::endl;
                    /*std::cout << "Interseções: " << intersection_points.size() << " " << std::endl;
                    print_points(intersection_points);
                    std::cout << "Normais de face: " << face_normals.size() << " " << std::endl;
                    print_points(face_normals);
                    std::cout << "Sinal P1: " << sign1 << " " << "Near P1: " << nearmost_index1 << " ";
                    print_point(p1);
                    std::cout << std::endl;
                    std::cout << "Sinal P2: " << sign2 << " " << "Near P2: " << nearmost_index2 << " ";
                    print_point(p2);
                    std::cout << std::endl;
                    std::cout << "Sinal Corretor: " << corrector << std::endl;*/
                    //TODO: investigate the reason of the divergence
                }
            }
            else
            {
                sign1 = computeSideOfPoint(p1, intersection_points[0], face_normals[0]);
                sign2 = sign1 == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;

                int corrector = computeSideOfPoint(p2, intersection_points[0], face_normals[0]);
                if (corrector != sign2){
                    std::cout << "DIVERGENCIA NO CALCULO DOS SINAIS EM CASO SINGLE" << std::endl;
                    /*std::cout << "Interseções: " << intersection_points.size() << " " << std::endl;
                    print_points(intersection_points);
                    std::cout << "Normais de face: " << face_normals.size() << " " << std::endl;
                    print_points(face_normals);
                    std::cout << "Sinal P1: " << sign1 << " ";
                    print_point(p1);
                    std::cout << std::endl;
                    std::cout << "Sinal P2: " << sign2 << " ";
                    print_point(p2);
                    std::cout << std::endl;
                    std::cout << "Sinal Corretor: " << corrector << std::endl;*/
                    //TODO: investigate the reason of the divergence
                }
            }

            if (vecsigns[c1] == MATERIAL_UNKNOWN){
                vecsigns[c1] = sign1;
            }
            else{
                if (vecsigns[c1] != sign1){
                    std::cout << intersection_points.size() << " edge1 " << c1 << "--" << c2 << " : " << vecsigns[c1] << " --> " << sign1 << std::endl;
                    inversion_flag = true;
                    //should_add = false;




                    //if (sign1 != vecsigns[c2])
                    vecsigns[c1] = sign1;

                }
            }
            if (vecsigns[c2] == MATERIAL_UNKNOWN){
                //vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
                vecsigns[c2] = sign2;

            }
            else{
                if (vecsigns[c2] != sign2){
                    std::cout << intersection_points.size() << " edge2 " << c2 << "--" << c1 << " : " << vecsigns[c2] << " --> " << sign2 << std::endl;
                    inversion_flag = true;
                    //should_add = false;

                    //if (sign2 != vecsigns[c1])
                    vecsigns[c2] = sign2;
                }
            }
            /*END VERTEX CLASSIFICATION*/
            edges_intersected++;
            /*std::cout << "Classificação das arestas. Aresta: " << c1 << "--" << c2 << std::endl;
            std::cout << "Interseções: " << intersection_points.size() << " " << std::endl;
            print_points(intersection_points);
            std::cout << "Normais de face: " << face_normals.size() << " " << std::endl;
            print_points(face_normals);
            std::cout << "Sinal P1: " << sign1 << " " << "Vecsigns[c1]: " << vecsigns[c1];
            print_point(p1);
            std::cout << std::endl;
            std::cout << "Sinal P2: " << sign2 << " " << "Vecsigns[c2]: " << vecsigns[c2];
            print_point(p2);
            std::cout << std::endl;*/

            std::string ehash = hashedge(p1, p2);
            HermiteData data;
            if (should_add) {
                //vec3 &n = normals[0];
                vec3 &n = face_normals[0];
                vec3 &v = intersection_points[0];
                qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
                averageNormal += n;
                hasIntersection = true;

                leaf->num_intersections++;
                data.intersection = intersection_points[0];
                //data.normal = normals[0];
                data.normal = face_normals[0];
            }
            edgepool[ehash] = data;
        }
    }

    if (inversion_flag){
        std::cout << "Current Signs: ";
        for (int j = 0; j < 8; ++j) {
            std::cout << vecsigns[j];
        }
        std::cout << std::endl;
        std::cout << "Old Signs: ";
        for (int j = 0; j < 8; ++j) {
            auto vertex = leaf->get_vertex(j);
            int sign = leafvertexpool[hashvertex(vertex)];
            std::cout << sign;
        }
        std::cout << std::endl;
        //exit(89);

        std::ofstream inversionfile;
        std::cout << "Num edges intersected: " << edges_intersected << std::endl;
        inversionfile.open("../subproducts/inversion.ply", std::ios::app);
        for (int i = 0; i < 8; ++i) {
            const vec3 cornerPos = leaf->get_vertex(i);
            int sign = vecsigns[i];
            if (sign == MATERIAL_SOLID) {
                inversionfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 255 << " " << 0<< std::endl;
            }

            if (sign == MATERIAL_AIR) {
                inversionfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 0 << " " << 255 << " " << 0 << std::endl;
            }

            if (sign == MATERIAL_UNKNOWN) {
                inversionfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
            }

        }
        inversionfile.close();
        //exit(62);
    }
    if (!hasIntersection){
        if (leaf->type == NODE_LEAF)
        {
            return leaf;
        }
        delete leaf;
        return nullptr;
    }

    int corners = 0;

    //updateSignsArray(vecsigns, 8);
    //updateSignsArray(vecsigns, 8, edges_intersected, leaf);
    //mergeSigns(vecsigns, leaf);

    /*for (size_t i = 0; i < 8; i++)
    {   //encode the signs to the corners variable to save memory
        if (vecsigns[i] == MATERIAL_AIR)
            corners |= (vecsigns[i] << i);
    }*/

    if (leaf->drawInfo == nullptr){
        leaf->drawInfo = new OctreeDrawInfo();
    }

    //leaf->drawInfo->corners |= corners;


    //TODO: pass drawinfo data from parent to child when descending
    //CHECK WHICH QEF TO KEEP

    bool old_opt = is_good_optimization(leaf, leaf->drawInfo->qef);
    bool cur_opt = is_good_optimization(leaf, qef.getData());
    if (old_opt && !cur_opt)
    {
        tmpQef.setData(leaf->drawInfo->qef);
    }
    else if (!old_opt && cur_opt){
        tmpQef.setData(qef.getData());
    }
    else {
        tmpQef.setData(leaf->drawInfo->qef + qef.getData());
    }



    //END CHECK

    //leaf->drawInfo->qef = leaf->drawInfo->qef + qef.getData();
    leaf->drawInfo->qef = tmpQef.getData();
    leaf->drawInfo->averageNormal += averageNormal;
    //leaf->drawInfo = drawInfo;
    leaf->type = NODE_LEAF;
    for (int i = 0; i < NUM_CHILDREN; ++i) {
        std::string vertex_hash = hashvertex(leaf->get_vertex(i));
        if (leafvertexpool.count(vertex_hash) == 0)
        {
            leafvertexpool[vertex_hash] = vecsigns[i];
        }
        else
        {
            int oldsign = leafvertexpool[vertex_hash];
            if (oldsign == MATERIAL_UNKNOWN){
                leafvertexpool[vertex_hash] = vecsigns[i];
                continue;
            }
            if (oldsign == MATERIAL_AMBIGUOUS){
                continue;
            }

            if (vecsigns[i] != MATERIAL_UNKNOWN && oldsign != MATERIAL_AMBIGUOUS){
                if (oldsign != vecsigns[i]){
#ifdef DEBUG
                    std::cout << " Computed Before: " << oldsign << " " << " Now: " << vecsigns[i] << std::endl;
                    Octree::divergence++;
#endif

                    // if we have divergence, we let the camera method classify
                    //leafvertexpool[vertex_hash] = MATERIAL_AMBIGUOUS;


                    //NOW WE CONSIDER OUT
                    leafvertexpool[vertex_hash] = vecsigns[i];//MATERIAL_AIR;
                }
            }
        }
    }

    //return clean_node(leaf);
    return leaf;
}

/*void Octree::classify_leaves_vertices(glm::vec3 cam_origin, OctreeNode* node, DefaultMesh &mesh)
{
    //trace("classify leaves vertices");
    if (node == nullptr) return;

    if (node->type == NODE_LEAF)
    {
        //trace("leaf");

        if (node->drawInfo == nullptr){
            node->drawInfo = new OctreeDrawInfo();
        }
        int corners = 0;
        for (int i = 0; i < NUM_CHILDREN; ++i)
        {
            vec3 cell_vertex = node->get_vertex(i);
            std::string vertex_hash = hashvertex(cell_vertex);
            if (leafvertexpool[vertex_hash] == MATERIAL_AMBIGUOUS)
            {
                //trace("updating pool");
                int sign = classify_vertex(cam_origin, cell_vertex, this->root, mesh);
                leafvertexpool[vertex_hash] = sign;
#ifdef DEBUG
                Octree::ambiguous_vertices++;
#endif

            }
            corners |= (leafvertexpool[vertex_hash] << i);
        }
        node->drawInfo->corners |= corners;
        //if (node->drawInfo->corners != corners){
        //   std::cout << "From intersection: " << node->drawInfo->corners << " From camera: " << corners << std::endl;
        //}
    }
    else
    {
        //trace("INTERNAL");
        for (int i = 0; i < NUM_CHILDREN; ++i) {
            classify_leaves_vertices(cam_origin, node->children[i], mesh);
        }
    }
}*/

// -------------------------------------------------------------------------------
void Octree::classify_leaves_vertices(OctreeNode* node)
{
    if (node == nullptr) return;

    if (node->type == NODE_INTERNAL)
    {
        for (int i = 0; i < NUM_CHILDREN; ++i)
        {
            classify_leaves_vertices(node->children[i]);
        }
    }

    if (node->type == NODE_LEAF)
    {
        //COMPUTING CORRECT VERTEX CLASSIFICATION
        int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                           MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
        bool should_revise_signs = false;
        for (int j = 0; j < 8; ++j)
        {
            vec3 vertex = node->get_vertex(j);
            std::string vertexhash = hashvertex(vertex);
            if (Octree::leafvertexpool.count(vertexhash) == 0){
                std::cout << "XABUZACO!!!!!!!! Contouring" << std::endl;
            }
            else{
                int stored_sign = Octree::leafvertexpool[vertexhash];

                if (stored_sign == MATERIAL_UNKNOWN || stored_sign == MATERIAL_AMBIGUOUS){
                    should_revise_signs = true;
                }
                else
                {
                    vecsigns[j] = stored_sign;
                }
            }
        }
        //TODO: check if i can save the sign in the hash during update to verify the ambiguity ratio
        if (should_revise_signs)
        {
            //std::cout << "pizza" << std::endl;
            updateSignsArray(vecsigns, 8, node);
            //updateSignsArray(vecsigns, 8);
            //std::cout << "hambuerguer" << std::endl;
        }
        for (int k = 0; k < 8; ++k) {
            node->drawInfo->corners |= (vecsigns[k] << k);
        }
        //COMPUTING CORRECT VERTEX CLASSIFICATION
    }
}

// -------------------------------------------------------------------------------
OctreeNode* Octree::SimplifyOctree(OctreeNode* node, const float threshold)
{
    if (!node)
    {
        //std::cout << "Empty node" << std::endl;
        return nullptr;
    }

    if (node->type != NODE_INTERNAL)
    {
        // can't simplify!
        return node;
    }
    //std::cout << "Simplifying at level: " << node->depth << std::endl;
    svd::QefSolver qef;

    int signs[8] = { MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                     MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN };
    int midsign = MATERIAL_UNKNOWN;
    int edgeCount = 0;
    bool isCollapsible = true;

    for (int i = 0; i < 8; ++i)
    {
        node->children[i] = SimplifyOctree(node->children[i], threshold);
        if (node->children[i])
        {
            OctreeNode* child = node->children[i];
            if (child->type == NODE_INTERNAL)
            {
                isCollapsible = false;
            }
            else
            {
                qef.add(child->drawInfo->qef);

                midsign = (child->drawInfo->corners >> (7 - i)) & 1;
                signs[i] = (child->drawInfo->corners >> i) & 1;

                edgeCount++;
            }
        }
    }

    if (!isCollapsible)
    {
        // at least one child is an internal node, can't collapse
        return node;
    }

    svd::Vec3 qefPosition;
    qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
    float error = qef.getError();

    // convert to glm vec3 for ease of use
    vec3 position(qefPosition.x, qefPosition.y, qefPosition.z);

    // at this point the masspoint will actually be a sum, so divide to make it the average
    if (error > threshold)
    {
        // this collapse breaches the threshold
        return node;
    }

    if ((position.x < node->min.x) || (position.x > (node->min.x + node->size)) ||
        (position.y < node->min.y) || (position.y > (node->min.y + node->size)) ||
        (position.z < node->min.z) || (position.z > (node->min.z + node->size)))
    {
        const auto& mp = qef.getMassPoint();
        position = vec3(mp.x, mp.y, mp.z);
    }

    // change the node from an internal node to a 'pseudo leaf' node
    OctreeDrawInfo* drawInfo = new OctreeDrawInfo;

    for (int i = 0; i < 8; i++)
    {
        if (signs[i] == MATERIAL_UNKNOWN)
        {
            // Undetermined, use centre sign instead
            drawInfo->corners |= (midsign << i);
        }
        else
        {
            drawInfo->corners |= (signs[i] << i);
        }
    }

    drawInfo->averageNormal = vec3(0.f);
    for (int i = 0; i < 8; ++i)
    {
        if (node->children[i])
        {
            OctreeNode* child = node->children[i];
            if (child->type != NODE_INTERNAL)
            {
                drawInfo->averageNormal += child->drawInfo->averageNormal;
            }
        }
    }

    drawInfo->averageNormal = glm::normalize(drawInfo->averageNormal);
    drawInfo->position = position;
    drawInfo->qef = qef.getData();

    for (int i = 0; i < 8; i++)
    {
        delete node->children[i];
        node->children[i] = nullptr;
    }

    node->type = NODE_PSEUDO;
    node->drawInfo = drawInfo;

    return node;
}
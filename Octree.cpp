//
// Created by Hallison da Paz on 18/11/2016.
//

#include <fstream>
#include "Octree.h"
#include "Utils.h"

#include "Reconstruction.h"

using glm::vec3;
// ----------------------------------------------------------------------------
std::unordered_map<std::string, int> Octree::leafvertexpool;
int Octree::unoptimized_points = 0;
int Octree::divergence = 0;
int Octree::ambiguous_vertices = 0;
// ----------------------------------------------------------------------------
int classify_vertex(glm::vec3 vertex, glm::vec3 cam_origin, OctreeNode* root, DefaultMesh &mesh);
// ----------------------------------------------------------------------------


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
// TODO: check if this can help
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
    trace("classifying vertices");

    classify_leaves_vertices(cam_origin, this->root, mesh);
    std::cout << "Divergence: " << Octree::divergence << std::endl;
    std::cout << "Ambiguities solved: " << Octree::ambiguous_vertices << std::endl;
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
    if (((node->parent && node->parent->innerEmpty()) || node->depth == max_depth) && node->depth > 4)
    {
        //return ConstructLeaf(node, max_depth, mesh);
        return ConstructLeafIntersection(node, max_depth, mesh);
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
        if (node->depth < max_depth && !node->parent->innerEmpty())
        {
            node->type = NODE_INTERNAL;
            node->construct_or_update_children(max_depth, mesh);
        }
        else
        {
            //node = update_leaf(node, max_depth, mesh);
            node = update_leaf_intersection(node, max_depth, mesh);
        }
        return node;
    }

    node->construct_or_update_children(max_depth, mesh);
    return node;
}

OctreeNode *Octree::ConstructLeaf(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh) {
    if (!leaf)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }

    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    bool hasIntersection = false;
    //TODO: optimize computations to avoid redundant intersections (same edge from other cell)
    for (int i = 0; i < 12; ++i) //for each edge
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

        vec3 intersection;
        std::vector<vec3> intersection_points, normals;//, face_normals;
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
                if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                    continue;
                }
                intersection_points.push_back(intersection);

                float u, v, w;
                barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);
                vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                normal_at_intersection =  glm::normalize(normal_at_intersection);
//                normals.push_back(normal_at_intersection);
                vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                //face_normals.push_back(face_normal);
                //hasIntersection = true;
                normals.push_back(face_normal);
            }
        }
        if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
            if (leaf->depth < max_depth){
                std::cout << intersection_points.size() << " Child Depth: " << leaf->depth+1 << " Child Size: " << leaf->size/2 << std:: endl;

                //leaf->type = NODE_INTERNAL;
                if(leaf->construct_or_update_children(max_depth, mesh))
                {
                    leaf->type = NODE_INTERNAL;
                    return leaf;
                }
                std::cout << "SERIAO????" << std::endl; //if it has an intersection why not the children?
                delete leaf;
                return nullptr;
            }
        }
        // if we consider that an intersection happened.
        // we'll consider only the first intersection for now
        if (intersection_points.size() > 0)
        {
            vec3 &n = normals[0];
            vec3 &v = intersection_points[0];
            qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
            averageNormal += n;
            hasIntersection = true;
        }
    }

    if (!hasIntersection)
    {   // voxel is full inside or outside the volume
        delete leaf;
        return nullptr;
    }

    if (leaf->drawInfo == nullptr){
        leaf->drawInfo = new OctreeDrawInfo();
    }
    leaf->drawInfo->qef = qef.getData();
    leaf->drawInfo->averageNormal += averageNormal;
    //leaf->drawInfo = drawInfo;
    leaf->type = NODE_LEAF;
    for (int i = 0; i < NUM_CHILDREN; ++i) {
        std::string vertex_hash = hashvertex(leaf->get_vertex(i));
        if (leafvertexpool.count(vertex_hash) == 0)
            leafvertexpool[vertex_hash] = MATERIAL_UNKNOWN;
    }
    //return clean_node(leaf);
    return leaf;
}

OctreeNode *Octree::update_leaf(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh)
{
    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    bool hasIntersection = false;
    //std::cout << "So far, so good" << std::endl;
    //TODO: optimize computations to avoid redundant intersections (same edge from other cell)
    for (int i = 0; i < NUM_EDGES; ++i) //for each edge
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

        vec3 intersection;
        std::vector<vec3> intersection_points, normals;//, face_normals;
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
                if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                    continue;
                }
                intersection_points.push_back(intersection);

                float u, v, w;
                barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);
                vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                normal_at_intersection =  glm::normalize(normal_at_intersection);
                normals.push_back(normal_at_intersection);
                //vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                //face_normals.push_back(face_normal);
                //hasIntersection = true;
            }
        }
        if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
            if (leaf->depth < max_depth){
                std::cout << intersection_points.size() << " Child Depth: " << leaf->depth+1 << " Child Size: " << leaf->size/2 << std:: endl;

                //leaf->type = NODE_INTERNAL;
                if(leaf->construct_or_update_children(max_depth, mesh))
                {
                    leaf->type = NODE_INTERNAL;
                    return leaf;
                }
                std::cout << "SERIAO????" << std::endl; //if it has an intersection why not the children?
                return leaf;
            }
        }
        // if we consider that an intersection happened.
        // we'll consider only the first intersection for now
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
    //return clean_node(leaf);
    return leaf;
}

OctreeNode *Octree::update_leaf_intersection(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh)
{
    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    bool hasIntersection = false;
    //std::cout << "So far, so good" << std::endl;
    int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                       MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
    //TODO: optimize computations to avoid redundant intersections (same edge from other cell)
    for (int i = 0; i < NUM_EDGES; ++i) //for each edge
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

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
                if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                    continue;
                }
                intersection_points.push_back(intersection);

                float u, v, w;
                barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);
                vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                normal_at_intersection =  glm::normalize(normal_at_intersection);
//                normals.push_back(normal_at_intersection);
                vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                face_normals.push_back(face_normal);
                //hasIntersection = true;
                normals.push_back(face_normal);
            }
        }
        if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
            if (leaf->depth < max_depth){
                std::cout << intersection_points.size() << " Child Depth: " << leaf->depth+1 << " Child Size: " << leaf->size/2 << std:: endl;

                //leaf->type = NODE_INTERNAL;
                if(leaf->construct_or_update_children(max_depth, mesh))
                {
                    leaf->type = NODE_INTERNAL;
                    return leaf;
                }
                std::cout << "SERIAO????" << std::endl; //if it has an intersection why not the children?
                return leaf;
            }
        }
        // if we consider that an intersection happened.
        // we'll consider only the first intersection for now
        if (intersection_points.size() > 0)
        {
            vec3 &n = normals[0];
            vec3 &v = intersection_points[0];
            qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
            averageNormal += n;
            hasIntersection = true;

            /*BEGIN VERTEX CLASSIFICATION*/
            int sign1 = computeSideOfPoint(p1, intersection_points[0], face_normals[0]);
            if (vecsigns[c1] == MATERIAL_UNKNOWN){
                vecsigns[c1] = sign1;
            }
            else{
                if (vecsigns[c1] != sign1){
                    std::cout << vecsigns[c1] << " --> " << sign1 << std::endl;
                }
            }
            if (vecsigns[c2] != MATERIAL_SOLID/*vecsigns[c2] == MATERIAL_UNKNOWN*/){
                vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
            }
            /*END VERTEX CLASSIFICATION*/
        }
    }

    if (!hasIntersection){
        return leaf;
    }

    int corners = 0;

    updateSignsArray(vecsigns, 8);
    //mergeSigns(vecsigns, leaf);

    for (size_t i = 0; i < 8; i++)
    {   //encode the signs to the corners variable to save memory
        corners |= (vecsigns[i] << i);
    }

    if (leaf->drawInfo == nullptr){
        leaf->drawInfo = new OctreeDrawInfo();
    }

    leaf->drawInfo->corners |= corners;
    //TODO: pass drawinfo data from parent to child when descending
    leaf->drawInfo->qef = leaf->drawInfo->qef + qef.getData();
    leaf->drawInfo->averageNormal += averageNormal;
    //leaf->drawInfo = drawInfo;
    leaf->type = NODE_LEAF;
    for (int i = 0; i < NUM_CHILDREN; ++i) {
        std::string vertex_hash = hashvertex(leaf->get_vertex(i));
        if (leafvertexpool.count(vertex_hash) == 0)
        {
            //leafvertexpool[vertex_hash] = MATERIAL_UNKNOWN;
            leafvertexpool[vertex_hash] = vecsigns[i];
        }
        else
        {
            int oldsign = leafvertexpool[vertex_hash];
            assert(oldsign != MATERIAL_UNKNOWN); //if we run updateSignsArray

            if (vecsigns[i] != MATERIAL_UNKNOWN && oldsign != MATERIAL_AMBIGUOUS){
                if (oldsign != vecsigns[i]){
                    std::cout << "Computed Before: " << oldsign << " " << " Now: " << vecsigns[i] << std::endl;
                    Octree::divergence++;
                    // if we have divergence, we let the camera method classify
                    leafvertexpool[vertex_hash] = MATERIAL_AMBIGUOUS;
                }
            }
        }
    }

    //return clean_node(leaf);
    return leaf;
}


OctreeNode *Octree::ConstructLeafIntersection(OctreeNode *leaf, unsigned int max_depth, const DefaultMesh &mesh) {
    if (!leaf)
    {
        std::cout << "Trying to construct a leaf in the middle" << std::endl;
        return nullptr;
    }

    // otherwise the voxel contains the surface, so find the edge intersections
    vec3 averageNormal(0.f);
    svd::QefSolver qef;
    bool hasIntersection = false;
    int vecsigns[8] = {MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN,
                       MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN, MATERIAL_UNKNOWN};
    //TODO: optimize computations to avoid redundant intersections (same edge from other cell)
    for (int i = 0; i < 12; ++i) //for each edge
    {
        const int c1 = edgevmap[i][0];
        const int c2 = edgevmap[i][1];
        const vec3 p1 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c1]);
        const vec3 p2 = vec3(leaf->min + leaf->size*CHILD_MIN_OFFSETS[c2]);

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
                if ((intersection_points.size() > 0) && (glm::distance(intersection, intersection_points[0]) < POINT_DISTANCE_THRESHOLD)){
                    continue;
                }
                intersection_points.push_back(intersection);

                float u, v, w;
                barycentric(intersection, face_vertices[0], face_vertices[1], face_vertices[2], u, v, w);
                vec3 normal_at_intersection = u * openmesh_to_glm(mesh.normal(a)) + v * openmesh_to_glm(mesh.normal(b)) + w * openmesh_to_glm(mesh.normal(c));
                normal_at_intersection =  glm::normalize(normal_at_intersection);
//                normals.push_back(normal_at_intersection);
                vec3 face_normal = openmesh_to_glm(mesh.normal(*face));
                face_normals.push_back(face_normal);
                //hasIntersection = true;
                normals.push_back(face_normal);
            }
        }
        if (intersection_points.size() > 1) {
//            std::cout << intersection_points.size() << " Interseções na mesma aresta " << vecsigns[c1] << vecsigns[c2] << std::endl;
            if (leaf->depth < max_depth){
                std::cout << intersection_points.size() << " Child Depth: " << leaf->depth+1 << " Child Size: " << leaf->size/2 << std:: endl;

                //leaf->type = NODE_INTERNAL;
                if(leaf->construct_or_update_children(max_depth, mesh))
                {
                    leaf->type = NODE_INTERNAL;
                    return leaf;
                }
                std::cout << "SERIAO????" << std::endl; //if it has an intersection why not the children?
                delete leaf;
                return nullptr;
            }
        }
        // if we consider that an intersection happened.
        // we'll consider only the first intersection for now
        if (intersection_points.size() > 0)
        {
            vec3 &n = normals[0];
            vec3 &v = intersection_points[0];
            qef.add(v.x, v.y, v.z, n.x, n.y, n.z);
            averageNormal += n;
            hasIntersection = true;

            /*BEGIN VERTEX CLASSIFICATION*/
            int sign1 = computeSideOfPoint(p1, intersection_points[0], face_normals[0]);
            if (vecsigns[c1] == MATERIAL_UNKNOWN){
                vecsigns[c1] = sign1;
            }
            else{
                if (vecsigns[c1] != sign1){
                    std::cout << vecsigns[c1] << " --> " << sign1 << std::endl;
                }
            }
            if (vecsigns[c2] != MATERIAL_SOLID/*vecsigns[c2] == MATERIAL_UNKNOWN*/){
                vecsigns[c2] = vecsigns[c1] == MATERIAL_AIR ? MATERIAL_SOLID : MATERIAL_AIR;
            }
            /*END VERTEX CLASSIFICATION*/
        }
    }

    if (!hasIntersection)
    {   // voxel is full inside or outside the volume
        delete leaf;
        return nullptr;
    }

    int corners = 0;
    updateSignsArray(vecsigns, 8);
    //mergeSigns(vecsigns, leaf);

    for (size_t i = 0; i < 8; i++)
    {   //encode the signs to the corners variable to save memory
        corners |= (vecsigns[i] << i);
        //updateVertexpool(OctreeNode::vertexpool, leaf->min + leaf->size*CHILD_MIN_OFFSETS[i], vecsigns[i]);
    }

    if (leaf->drawInfo == nullptr){
        leaf->drawInfo = new OctreeDrawInfo();
    }

    leaf->drawInfo->corners |= corners;

    leaf->drawInfo->qef = qef.getData();
    leaf->drawInfo->averageNormal += averageNormal;
    //leaf->drawInfo = drawInfo;
    leaf->type = NODE_LEAF;
    for (int i = 0; i < NUM_CHILDREN; ++i) {
        std::string vertex_hash = hashvertex(leaf->get_vertex(i));
        if (leafvertexpool.count(vertex_hash) == 0)
        {
            //leafvertexpool[vertex_hash] = MATERIAL_UNKNOWN;
            leafvertexpool[vertex_hash] = vecsigns[i];
        }
        else
        {
            int oldsign = leafvertexpool[vertex_hash];
            assert(oldsign != MATERIAL_UNKNOWN); //if we run updateSignsArray
            /*if (oldsign == MATERIAL_UNKNOWN){
                leafvertexpool[vertex_hash] = vecsigns[i];
                std::cout << "XABUZACO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            }
            if (oldsign == MATERIAL_AMBIGUOUS){
                std::cout << "AMBIGUITY HERE!!!!!" << std::endl;
                continue;
            }*/
            if (vecsigns[i] != MATERIAL_UNKNOWN && oldsign != MATERIAL_AMBIGUOUS){
                if (oldsign != vecsigns[i]){
                    std::cout << "Computed Before: " << oldsign << " " << " Now: " << vecsigns[i] << std::endl;
                    Octree::divergence++;
                    // if we have divergence, we let the camera method classify
                    leafvertexpool[vertex_hash] = MATERIAL_AMBIGUOUS;
                }
            }
        }

    }
    //return clean_node(leaf);
    return leaf;
}

void Octree::classify_leaves_vertices(glm::vec3 cam_origin, OctreeNode* node, DefaultMesh &mesh)
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
            if (leafvertexpool[vertex_hash] == MATERIAL_AMBIGUOUS/*MATERIAL_UNKNOWN*/)
            {
                //trace("updating pool");
                int sign = classify_vertex(cam_origin, cell_vertex, this->root, mesh);
                leafvertexpool[vertex_hash] = sign;
                Octree::ambiguous_vertices++;

            }
            corners |= (leafvertexpool[vertex_hash] << i);
        }
        node->drawInfo->corners |= corners;
        /*if (node->drawInfo->corners != corners){
            std::cout << "From intersection: " << node->drawInfo->corners << " From camera: " << corners << std::endl;
        }*/
    }
    else
    {
        //trace("INTERNAL");
        for (int i = 0; i < NUM_CHILDREN; ++i) {
            classify_leaves_vertices(cam_origin, node->children[i], mesh);
        }
    }
}

// Intersect ray R(t) = p + t*d against AABB a. When intersecting,
// return intersection distance tmin and point q of intersection
bool intersectRayBox(vec3 origin, vec3 dest, const glm::vec3 min, const Real size, Real &intersection_distance, vec3 &intersection) {
    intersection_distance = 0.0f; // set to -FLT_MAX to get first hit on line
    float tmax = FLT_MAX; // set to max distance ray can travel (for segment)
    // For all three slabs
    vec3 dir = glm::normalize(dest - origin);
    vec3 max = min + vec3(size);
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(dir[i]) < EPSILON)
        {
            // Ray is parallel to slab. No hit if origin not within slab
            if (origin[i] < min[i] || origin[i] > max[i])
                return 0;
        }
        else
        {
            // Compute intersection t value of ray with near and far plane of slab
            float ood = 1.0f / dir[i];
            float t1 = (min[i] - origin[i]) * ood;
            float t2 = (max[i] - origin[i]) * ood;
            // Make t1 be intersection with near plane, t2 with far plane
            if (t1 > t2)
            {
                Real aux = t1;
                t1 = t2;
                t2 = aux;
            }
            // Compute the intersection of slab intersection intervals
            if (t1 > intersection_distance) intersection_distance = t1;
            if (t2 > tmax) tmax = t2;
            // Exit with no collision as soon as slab intersection becomes empty
            if (intersection_distance > tmax)
                return false;
        }
    }
    // Ray intersects all 3 slabs. Return point (q) and intersection t value (intersection_distance)
    intersection = origin + dir * intersection_distance;
    return true;
}


int ray_faces_intersection(const glm::vec3 origin, const glm::vec3 dest, DefaultMesh &mesh,
                           std::list<DefaultMesh::FaceHandle> &facelist, std::unordered_map<int, bool> &visited_triangles)
{
    int num_intersections = 0;
    std::vector<vec3> intersections;
    for (std::list<DefaultMesh::FaceHandle>::iterator face = facelist.begin(); face != facelist.end(); ++face){
        if (visited_triangles.count(face->idx()) == 0){
            auto fv_it = mesh.cfv_iter(*face);
            DefaultMesh::VertexHandle a = *fv_it;
            DefaultMesh::VertexHandle b = *(++fv_it);
            DefaultMesh::VertexHandle c = *(++fv_it);

            vec3 face_vertices[3] = {openmesh_to_glm(mesh.point(a)), openmesh_to_glm(mesh.point(b)), openmesh_to_glm(mesh.point(c))};
            Vertex vertices[3] = { face_vertices[0], face_vertices[1], face_vertices[2]};

            vec3 intersection;

            if (moller_triangle_intersection(origin, dest, vertices, intersection)) {
                //keeps the intersection here
                intersections.push_back(intersection);
                ++num_intersections;
            }
            visited_triangles[face->idx()] = true;
        }
    }
    /*if (num_intersections > 1){
        for (int i = 0; i < intersections.size(); ++i) {
            std::cout << intersections[i].x << " " <<  intersections[i].y << " " <<  intersections[i].z << std::endl;
        }
        exit(67);
    }*/
    return num_intersections;
}

int ray_mesh_intersection(glm::vec3 cam_origin, glm::vec3 vertex, OctreeNode* root, DefaultMesh &mesh, std::unordered_map<int, bool> &visited_triangles)
{
    if (root == nullptr){
        return 0;
    }
    int num_intersections = 0;
    vec3 intersection;
    Real t;
    if (intersectRayBox(cam_origin, vertex, root->min, root->size, t, intersection)){
        if (root->type == NODE_LEAF)
        {
            num_intersections += ray_faces_intersection(cam_origin, vertex, mesh, root/*->meshInfo*/->crossingFaces, visited_triangles);
            num_intersections += ray_faces_intersection(cam_origin, vertex, mesh, root/*->meshInfo*/->innerFaces, visited_triangles);
        }
        else
        {
            for (int i = 0; i < NUM_CHILDREN; ++i)
            {
                num_intersections += ray_mesh_intersection(cam_origin, vertex, root->children[i], mesh, visited_triangles);
            }
        }
    }
    return num_intersections;
}

int classify_vertex(glm::vec3 cam_origin, glm::vec3 vertex, OctreeNode* root, DefaultMesh &mesh)
{
    std::unordered_map<int, bool> visited_triangles;
    int num_intersections = ray_mesh_intersection(cam_origin, vertex, root, mesh, visited_triangles);
    if (num_intersections >= 2){
        //std::cout << "intersections: " << num_intersections << std::endl;
    }
    if (/*num_intersections%2 == 1*/num_intersections > 0)
    {
        return MATERIAL_SOLID;
    }
    return MATERIAL_AIR;
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
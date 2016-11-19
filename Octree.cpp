//
// Created by Hallison da Paz on 18/11/2016.
//

#include <fstream>
#include "Octree.h"
#include "Utils.h"


using glm::vec3;
// ----------------------------------------------------------------------------
std::unordered_map<std::string, int> Octree::leafvertexpool;
// ----------------------------------------------------------------------------
int classify_vertex(glm::vec3 vertex, glm::vec3 cam_origin, OctreeNode* root, DefaultMesh &mesh);
// ----------------------------------------------------------------------------


bool OctreeNode::construct_children(unsigned int max_depth, const DefaultMesh &mesh)
{
    const float childSize = this->size / 2;
    const int childHeight = this->depth + 1;
    bool hasChildren = false;
    for (int i = 0; i < NUM_CHILDREN; i++)
    {
        OctreeNode* child = new OctreeNode(NODE_INTERNAL,
                                           this->min + (CHILD_MIN_OFFSETS[i] * childSize),
                                           childSize,
                                           childHeight,
                                           this);
        this->children[i] = Octree::BuildMeshHierarchy(child, max_depth, mesh);
        hasChildren |= (this->children[i] != nullptr);
    }
    return hasChildren;
}

Octree::Octree(glm::vec3 min, Real size, unsigned int max_depth, DefaultMesh &mesh)
{
    root = new OctreeNode(NODE_INTERNAL, min, size, 0);
    BuildMeshHierarchy(root, max_depth, mesh);
    vec3 cam_origin(0, -8, 0);
    classify_leaves_vertices(cam_origin, this->root, mesh);
}

OctreeNode *Octree::BuildMeshHierarchy(OctreeNode *node, unsigned int max_depth, const DefaultMesh &mesh)
{
    if (!node) return nullptr;
    //trace("begin");
    /*if (node->meshInfo == nullptr)
        node->meshInfo = new OctreeMeshInfo();*/
    select_inner_crossing_faces(node, mesh);
    //trace("hop");
    if (node->innerEmpty() && node->crossingEmpty())
    {   //Empty space, no triangles crossing or inside this cell
        delete node;
        return nullptr;
    }
    if ((node->parent && node->parent->innerEmpty()) || node->depth == max_depth)
    {
        node->type = NODE_LEAF;
        for (int i = 0; i < NUM_CHILDREN; ++i) {
            std::string vertex_hash = hashvertex(node->get_vertex(i));
            if (leafvertexpool.count(vertex_hash) == 0)
                leafvertexpool[vertex_hash] = MATERIAL_UNKNOWN;
        }
        return node;
    }
    if (node->construct_children(max_depth, mesh))
    {
        return node;
    }

    delete node;
    return nullptr;
}

void Octree::classify_leaves_vertices(glm::vec3 cam_origin, OctreeNode* node, DefaultMesh &mesh)
{
    //trace("classify leaves vertices");
    if (node == nullptr) return;

    if (node->type == NODE_LEAF)
    {
        //trace("leaf");
        std::ofstream interiorfile, exteriorfile;
        interiorfile.open("../subproducts/interior_color_updated.ply", std::ios::app);
        exteriorfile.open("../subproducts/exterior_color_updated.ply", std::ios::app);
        node->drawInfo = new OctreeDrawInfo();
        for (int i = 0; i < NUM_CHILDREN; ++i)
        {
            vec3 cell_vertex = node->get_vertex(i);
            std::string vertex_hash = hashvertex(cell_vertex);
            if (leafvertexpool[vertex_hash] == MATERIAL_UNKNOWN)
            {
                //trace("updating pool");
                int sign = classify_vertex(cam_origin, cell_vertex, this->root, mesh);
                leafvertexpool[vertex_hash] = sign;
                // DEBUG --------------------------------------------------------------------
                const vec3 cornerPos = node->min + CHILD_MIN_OFFSETS[i]*node->size;
                if (sign == MATERIAL_SOLID) {
                    interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 128 << " " << 255 << std::endl;
                }

                if (sign == MATERIAL_AIR) {
                    exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 64 << " " << 255 << " " << 64 << std::endl;
                }

                if (sign == MATERIAL_UNKNOWN) {
                    interiorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
                    exteriorfile << cornerPos.x << " " << cornerPos.y << " " << cornerPos.z << " " << 255 << " " << 0 << " " << 0 << std::endl;
                }
            }
            node->drawInfo->corners |= (leafvertexpool[vertex_hash] << i);
        }
        interiorfile.close();
        exteriorfile.close();
        // DEBUG --------------------------------------------------------------------//
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
    if (num_intersections > 2){
        std::cout << "intersections: " << num_intersections << std::endl;
    }
    if (num_intersections%2 == 1)
    {
        return MATERIAL_SOLID;
    }
    return MATERIAL_AIR;
}
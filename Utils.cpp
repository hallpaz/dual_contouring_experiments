#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <list>

#include "Utils.h"
#include "Constants.h"
#include "glm/ext.hpp"

#include "old/density.h"
#include "Octree.h"

using glm::vec3;


// ----------------------------------------------------------------------------

void divideFacesByLocation(OctreeNode *node, std::list<DefaultMesh::FaceHandle> &facesList, const DefaultMesh &mesh)
{
    //trace("divide face by location");
    auto f_it = facesList.begin();
    //parent's inner triangles
    while(f_it != facesList.end())
    {
        if (! mesh.is_valid_handle(*f_it))
        {
            std::cout << "Invalid Face Handle on divideFacesByLocation" << std::endl;
            ++f_it;
            continue;
        }
        switch (triangleRelativePosition(mesh, *f_it, node->min, node->size))
        {
            case INSIDE:
                node/*->meshInfo*/->innerFaces.push_back(*f_it);
                /*If the triangle is located inside the cell, we remove it from the cell's parent list*/
                f_it = facesList.erase(f_it);
                break;
            case CROSSING:
                node/*->meshInfo*/->crossingFaces.push_back(*f_it);
                /*If the triangle might cross the cell, we can't remove it from the cell's parent list
                 * because it might cross other cells as well*/
                ++f_it;
                break;
            default:
                ++f_it;
        }
    }
}

// ----------------------------------------------------------------------------

void select_inner_crossing_faces(OctreeNode *node, const DefaultMesh &mesh)
{
    //trace("select inner crossing");
    node->clean();
    if (node->parent != nullptr)
    {
        divideFacesByLocation(node, node->parent/*->meshInfo*/->innerFaces, mesh); //parent's inner triangles
        divideFacesByLocation(node, node->parent/*->meshInfo*/->crossingFaces, mesh); //parent's crossing triangles
    }
    else
    {   //initializes the parent list with all triangles
        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
        {
            node/*->meshInfo*/->innerFaces.push_back(*f_it);
        }
    }
}

// ----------------------------------------------------------------------------

int computeSideOfPoint(const glm::vec3 point, const glm::vec3 intersection, const glm::vec3 face_normal)
{
    return glm::dot(point - intersection, face_normal) < 0.f ? MATERIAL_SOLID : MATERIAL_AIR;
}
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------
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

// -------------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
void updateVertexpool(std::unordered_map<std::string, int> &pool, const glm::vec3 &vertex, int &sign)
{
    if(sign == MATERIAL_UNKNOWN)
    {
        std::cout << "TENTA OUTRA NEGAO!!" << std::endl;
        exit(98);
    }
    std::string vertexhash = hashvertex(vertex);
    if (pool.count(vertexhash) == 0)
    {
        pool[vertexhash] = sign;
    }
    else
    {
        if (pool[vertexhash] == MATERIAL_SOLID)
        {   // if it was marked as interior anytime
            sign = pool[vertexhash];
        }
        else
        {
            pool[vertexhash] = sign;
        }
    }
    if (sign != MATERIAL_UNKNOWN && pool[vertexhash] == MATERIAL_UNKNOWN){
        std::cout << "LEFTING UNKNOWN!!!!!!!" << std::endl;
        exit(98);
    }
}

// ----------------------------------------------------------------------------
void mergeSigns(int *vecsigns, OctreeNode* node)
{
    for (int i = 0; i < 8; ++i) {
        if (vecsigns[i] == MATERIAL_UNKNOWN){
            vec3 vertex = node->get_vertex(i);
            std::string vertexhash = hashvertex(vertex);
            if (Octree::leafvertexpool.count(vertexhash) != 0){
                int storedsign = Octree::leafvertexpool[vertexhash];
                if (storedsign != MATERIAL_AMBIGUOUS){
                    vecsigns[i] = storedsign;
                }
                else {
                    vecsigns[i] = 0;
                }
            }
            else{
                vecsigns[i] = 0;
            }
        }

    }
}

// ----------------------------------------------------------------------------
void updateSignsArray(int *vecsigns, int size)
{
    /*for (int k = 0; k < 8; ++k) {
        std::cout << vecsigns[k];
    }*/
    bool checksigns = true;
    while(checksigns) {
        checksigns = false;
        for (size_t i = 0; i < size; ++i) {
            if (vecsigns[i] == MATERIAL_UNKNOWN) {
                checksigns = true;
                for (int j = 0; j < 3; ++j) {
                    //DEBUG
                    /*bool flag = false;
                    if (vecsigns[vneighbors[i][0]] != vecsigns[vneighbors[i][1]]
                        && (vecsigns[vneighbors[i][0]] != MATERIAL_UNKNOWN && vecsigns[vneighbors[i][1]] != MATERIAL_UNKNOWN)){
                        flag = true;
                    }
                    if (vecsigns[vneighbors[i][0]] != vecsigns[vneighbors[i][2]]
                        && (vecsigns[vneighbors[i][0]] != MATERIAL_UNKNOWN && vecsigns[vneighbors[i][2]] != MATERIAL_UNKNOWN)){
                        flag = true;
                    }
                    if (vecsigns[vneighbors[i][1]] != vecsigns[vneighbors[i][2]]
                        && (vecsigns[vneighbors[i][1]] != MATERIAL_UNKNOWN && vecsigns[vneighbors[i][2]] != MATERIAL_UNKNOWN)){
                        flag = true;
                    }
                    if (flag) {
                        std::cout << i << ": Diferentes: " << vneighbors[i][0] << ": " << vecsigns[vneighbors[i][0]] << " "
                                << vneighbors[i][1] << ": " << vecsigns[vneighbors[i][1]] << " "
                                << vneighbors[i][2] << ": " << vecsigns[vneighbors[i][2]] << std::endl;
                        for (int k = 0; k < 8; ++k) {
                            std::cout << vecsigns[k];
                        }
                        std::cout << std::endl;
                    }*/
                    //DEBUG
                    int n = vneighbors[i][j];
                    if (vecsigns[n] != MATERIAL_UNKNOWN) {
                        vecsigns[i] = vecsigns[n];
                        break;
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
/*void updateSignsArray(int *vecsigns, int size, int edges_intersected, OctreeNode* node)
{

    for (int i = 0; i < NUM_CHILDREN; ++i) {
        if (vecsigns[i] == MATERIAL_UNKNOWN) {
            //checksigns = true;
            for (int j = 0; j < 3; ++j) {
                //DEBUG
                bool flag = false;
                if (vecsigns[vneighbors[i][0]] != vecsigns[vneighbors[i][1]]
                    && (vecsigns[vneighbors[i][0]] != MATERIAL_UNKNOWN &&
                        vecsigns[vneighbors[i][1]] != MATERIAL_UNKNOWN)) {
                    flag = true;
                }
                if (vecsigns[vneighbors[i][0]] != vecsigns[vneighbors[i][2]]
                    && (vecsigns[vneighbors[i][0]] != MATERIAL_UNKNOWN &&
                        vecsigns[vneighbors[i][2]] != MATERIAL_UNKNOWN)) {
                    flag = true;
                }
                if (vecsigns[vneighbors[i][1]] != vecsigns[vneighbors[i][2]]
                    && (vecsigns[vneighbors[i][1]] != MATERIAL_UNKNOWN &&
                        vecsigns[vneighbors[i][2]] != MATERIAL_UNKNOWN)) {
                    flag = true;
                }
                if (flag) {
                    std::cout << i << ": Diferentes: " << vneighbors[i][0] << ": " << vecsigns[vneighbors[i][0]] << " "
                              << vneighbors[i][1] << ": " << vecsigns[vneighbors[i][1]] << " "
                              << vneighbors[i][2] << ": " << vecsigns[vneighbors[i][2]] << std::endl;
                    for (int k = 0; k < 8; ++k) {
                        std::cout << vecsigns[k];
                    }
                    std::cout << " ei:" << edges_intersected << std::endl;
                    if (edges_intersected > 2) {
                        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << edges_intersected
                                  << std::endl;
                    }
                    node->irregular = true;
                    break;
                }
                //DEBUG
                int n = vneighbors[i][j];
                if (vecsigns[n] != MATERIAL_UNKNOWN) {
                    vecsigns[i] = vecsigns[n];
                    break;
                }
            }
        }

        std::string vertex_hash = hashvertex(node->get_vertex(i));
        if (Octree::leafvertexpool.count(vertex_hash) == 0)
        {
            Octree::leafvertexpool[vertex_hash] = vecsigns[i];
        }
        else
        {
            int stored_sign = Octree::leafvertexpool[vertex_hash];
            if (stored_sign == MATERIAL_UNKNOWN)
            {
                Octree::leafvertexpool[vertex_hash] = vecsigns[i];
            }
            else
            {
                if (vecsigns[i] != MATERIAL_UNKNOWN && stored_sign != MATERIAL_AMBIGUOUS){
                    if (stored_sign != vecsigns[i]){
                        std::cout << "Computed Before: " << stored_sign << " " << " Now: " << vecsigns[i] << std::endl;
                        Octree::divergence++;
                        // if we have divergence, we let the camera method classify
                        Octree::leafvertexpool[vertex_hash] = MATERIAL_AMBIGUOUS;
                    }
                }
            }
        }
    }*/
void updateSignsArray(int *vecsigns, int size, int edges_intersected, OctreeNode* node){
    bool checksigns = true;
    while(checksigns) {
        checksigns = false;
        for (size_t i = 0; i < size; ++i) {
            if (vecsigns[i] == MATERIAL_UNKNOWN) {
                checksigns = true;
                for (int j = 0; j < 3; ++j) {
                    //DEBUG
                    bool flag = false;
                    if (vecsigns[vneighbors[i][0]] != vecsigns[vneighbors[i][1]]
                        && (vecsigns[vneighbors[i][0]] != MATERIAL_UNKNOWN &&
                            vecsigns[vneighbors[i][1]] != MATERIAL_UNKNOWN)) {
                        flag = true;
                    }
                    if (vecsigns[vneighbors[i][0]] != vecsigns[vneighbors[i][2]]
                        && (vecsigns[vneighbors[i][0]] != MATERIAL_UNKNOWN &&
                            vecsigns[vneighbors[i][2]] != MATERIAL_UNKNOWN)) {
                        flag = true;
                    }
                    if (vecsigns[vneighbors[i][1]] != vecsigns[vneighbors[i][2]]
                        && (vecsigns[vneighbors[i][1]] != MATERIAL_UNKNOWN &&
                            vecsigns[vneighbors[i][2]] != MATERIAL_UNKNOWN)) {
                        flag = true;
                    }
                    if (flag) {
                        std::cout << i << ": Diferentes: " << vneighbors[i][0] << ": " << vecsigns[vneighbors[i][0]]
                                  << " "
                                  << vneighbors[i][1] << ": " << vecsigns[vneighbors[i][1]] << " "
                                  << vneighbors[i][2] << ": " << vecsigns[vneighbors[i][2]] << std::endl;
                        for (int k = 0; k < 8; ++k) {
                            std::cout << vecsigns[k];
                        }
                        std::cout << " ei:" << edges_intersected << std::endl;
                        if (edges_intersected > 2) {
                            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << edges_intersected
                                      << std::endl;
                        }
                        node->irregular = true;
                    }
                    //DEBUG
                    int n = vneighbors[i][j];
                    if (vecsigns[n] != MATERIAL_UNKNOWN) {
                        vecsigns[i] = vecsigns[n];
                        break;
                    }
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
 void updateSignsArray(int *vecsigns, int size, OctreeNode* node){
    bool checksigns = true;
//    for (int j = 0; j < 8; ++j) {
//        std::cout << vecsigns[j];
//    }
    bool terror = false;
    std::cout << std::endl;
    int index = 0;
    while(checksigns) {
        checksigns = false;
        index++;
        for (size_t i = 0; i < size; ++i) {
            if (vecsigns[i] == MATERIAL_UNKNOWN) {
                checksigns = true;
                int sn0 = vecsigns[vneighbors[i][0]];
                int sn1 = vecsigns[vneighbors[i][1]];
                int sn2 = vecsigns[vneighbors[i][2]];
                //std::cout << i << sn0 << sn1 << sn2 << std::endl;
                if (sn0 == MATERIAL_UNKNOWN && sn1 == MATERIAL_UNKNOWN && sn2 == MATERIAL_UNKNOWN){
                    std::cout << "caso do terror" << std::endl;
                    for (int j = 0; j < 8; ++j) {
                        std::cout << vecsigns[j];
                    }
                    std::cout << std::endl;
                    terror = true;
                    continue;
                }
                bool flag = false;
                if (sn0 != sn1 && (sn0 != MATERIAL_UNKNOWN && sn1 != MATERIAL_UNKNOWN)) {
                    flag = true;
                }
                if (sn0 != sn2 && (sn0 != MATERIAL_UNKNOWN && sn2 != MATERIAL_UNKNOWN)) {
                    flag = true;
                }
                if (sn1 != sn2 && (sn1 != MATERIAL_UNKNOWN && sn2 != MATERIAL_UNKNOWN)) {
                    flag = true;
                }
                if (flag) {
                    node->irregular = true;
                    vecsigns[i] = MATERIAL_AIR;
                }
                else
                {
                    vecsigns[i] = sn0 != MATERIAL_UNKNOWN ? sn0 : sn1 != MATERIAL_UNKNOWN ? sn1 : sn2;
                }
            }
        }
        if(index > 8){
            terror = true;
            break;
        }
    }
    if(terror){
        for (int j = 0; j < 8; ++j) {
            std::cout << vecsigns[j];
        }
        std::cout << std::endl;
        if(index > 8)
            exit(88);
    }
 }


// ----------------------------------------------------------------------------

void write_Ply(std::string filename, VertexBuffer &vertices, IndexBuffer &faces)
{
    std::ofstream outfile;
    outfile.open(filename.c_str(), std::ios::out);


    if(outfile.is_open()){
        outfile << "ply" << std::endl
                << "format ascii 1.0" << std::endl
                << "element vertex "<< vertices.size() << std::endl
                << "property float x" << std::endl
                << "property float y" << std::endl
                << "property float z" << std::endl
                << "property uchar red" << std::endl
                << "property uchar green" << std::endl
                << "property uchar blue" << std::endl
                << "element face " << faces.size() << std::endl
                << "property list uchar int vertex_index" << std::endl
                << "end_header" << std::endl;
        for (VertexBuffer::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) {
            outfile << v_it->position.x << " " << v_it->position.y << " " << v_it->position.z << " "
                    << v_it->color.x << " " << v_it->color.y << " " << v_it->color.z << std::endl;
        }

        for (IndexBuffer::iterator f_it = faces.begin(); f_it != faces.end(); ++f_it) {
            outfile << "3" << " " << f_it->a << " " << f_it->b << " " << f_it->c << std::endl;
        }

        outfile.flush();
        outfile.close();
    }
}


void write_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces)
{
    std::ofstream outfile;
    outfile.open(filename.c_str(), std::ios::out);

    if(outfile.is_open()){
        std::cout << "Arquivo aberto" << std::endl;

        outfile << "OFF" << std::endl;
        outfile << vertices.size() << " " << faces.size() << " 0" << std::endl;

        for ( std::vector<Vertex>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) {
            outfile << v_it->position.x << " " << v_it->position.y << " " << v_it->position.z << std::endl;
        }

        for (std::vector<Triangle>::iterator f_it = faces.begin(); f_it != faces.end(); ++f_it ) {
            outfile << "3" << " " << f_it->a << " " << f_it->b << " " << f_it->c << std::endl;
        }

        /*for (int i = 0; i < faces.size(); i+=3) {
            outfile << "3" << " " << faces[i] << " " << faces[i+1] << " " << faces[i+2] << std::endl;
        }*/

        outfile.flush();
        outfile.close();
    }
    else {
        std::cout << "Arquivo Fechado: ERRO" << std::endl;
    }

}

float read_OFF(std::string filename, std::vector<Vertex> &vertices, std::vector<Triangle> &faces, vec3 &minPoint) {
    std::ifstream inputfile;
    inputfile.open(filename.c_str(), std::ios::out);
    const float BIG = 999999.9;
    float minx = BIG, miny = BIG, minz = BIG;
    float maxx = -BIG, maxy = -BIG, maxz = -BIG;
    if (inputfile.is_open()){
        std::string off;
        inputfile >> off;
        if (off != "OFF"){
            std::cout << "Invalid OFF File" << std::endl;
        }
        int numVertices, numFaces, numEdges;
        inputfile >> numVertices >> numFaces >> numEdges;
        for (int i = 0; i < numVertices; ++i) {
            float x, y, z;
            inputfile >> x >> y >> z;
            if (x < minx) minx = x;
            if (y < miny) miny = y;
            if (z < minz) minz = z;
            if (x > maxx) maxx = x;
            if (y > maxy) maxy = y;
            if (z > maxz) maxz = z;
            vertices.push_back(Vertex(glm::vec3(x, y, z)));
        }
        for (int i = 0; i < numFaces; ++i) {
            int n, a, b, c;
            inputfile >> n >> a >> b >> c;
            faces.push_back(Triangle{a, b, c});
        }
        inputfile.close();
    }
    minPoint.x = minx;
    minPoint.y = miny;
    minPoint.z = minz;
    float xsize = maxx - minx;
    float ysize = maxy - miny;
    float zsize = maxz - minz;
    float size = xsize > ysize ? xsize : ysize;
    size = size > zsize ? size : zsize;

    //float size = fabs(maxd) > fabs(mind) ? ceilf(fabs(maxd)) : ceilf(fabs(mind));
    std::cout << "Bounding box Min: (" << minPoint.x << " " << minPoint.y << " " << minPoint.z << ") Size: " << size << std::endl;
    return size;
}

std::string hashvertex(const vec3& vertex){
    int precision = 3;
    std::stringstream rep;
    rep << vertex.x << std::setprecision (precision) << std::fixed
    << vertex.y << std::setprecision (precision) << std::fixed
    << vertex.z << std::setprecision (precision) << std::fixed;
    return rep.str();
}

std::string hashedge(const vec3& v1, const vec3& v2){
    float threshold = 0.0000001;
    if (fabs(v1.x - v2.x) > threshold){
        if (v1.x < v2.x){
            return (hashvertex(v1) + hashvertex(v2));
        }
        return (hashvertex(v2) + hashvertex(v1));
    }
    if (fabs(v1.y - v2.y) > threshold){
        if (v1.y < v2.y){
            return (hashvertex(v1) + hashvertex(v2));
        }
        return (hashvertex(v2) + hashvertex(v1));
    }

    if (v1.z < v2.z){
        return (hashvertex(v1) + hashvertex(v2));
    }
    return (hashvertex(v2) + hashvertex(v1));
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
// Implementation from the book Real Time Collision Detection
void barycentric(vec3 p, vec3 a, vec3 b, vec3 c, float &u, float &v, float &w)
{
    vec3 v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
    if (u < 0) u = 0.0;
}

bool moller_triangle_intersection(vec3 v1, vec3 v2, Vertex* triangle_vertices, vec3& intersection_point)
{
    float EPSILON = 0.000001;

    vec3 e1 = triangle_vertices[1].position - triangle_vertices[0].position;
    vec3 e2 = triangle_vertices[2].position - triangle_vertices[0].position;

    vec3 D = v2 - v1;
    auto P = glm::cross(D, e2);

    float det = glm::dot(e1, P);
    if (glm::abs(det) < EPSILON){
        return false;
    }

    float inv_det = 1.0/det;
    vec3 T = v1 - triangle_vertices[0].position;
    float u = glm::dot(T, P) * inv_det;

    if(u < 0.0f or u > 1.0f){
        return false;
    }

    vec3 Q = glm::cross(T, e1);
    float v = glm::dot(D, Q) * inv_det;
    if(v < 0.0f or (u + v)  > 1.0f) {
        return false;
    }

    float t = glm::dot(e2, Q) * inv_det;

    if(t > 0.0 && t < 1.0){
        intersection_point = v1 + (v2-v1)*t;
        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
/* This function was modified to use a binary search like solution.
 * */
vec3 ApproximateZeroCrossingPosition(const vec3& p0, const vec3& p1)
{
    const float tolerance = 0.00001;
    const int maxIterations = 21;
    float tbegin = 0.0;
    float tend = 1.0;
    float t = (tbegin + tend)/2;
    int currentIt = 0;
    float middle = Density_Func(p0 + ((p1 - p0) * t));
    float begin = Density_Func(p0 + ((p1 - p0) * tbegin));
    float end = Density_Func(p0 + ((p1 - p0) * tend));
    while (currentIt < maxIterations)
    {
        if (fabs(middle) < fabs(tolerance)){
            break;
        }
        //double check
        if (begin * end > 0){
            std::cout << "Weird, but this is not a bipolar edge! it: " << currentIt << std::endl;
        }
        //check the side of the middle point
        if(middle*begin > 0){ //same side of begin
            tbegin = t;
        }
        else { //same side of end
            tend = t;
        }
        t = (tbegin + tend)/2;
        middle = Density_Func(p0 + ((p1 - p0) * t));
        begin = Density_Func(p0 + ((p1 - p0) * tbegin));
        end = Density_Func(p0 + ((p1 - p0) * tend));
    }
    std::cout << "Approximation: " << middle << std::endl;
    return (p0 + ((p1 - p0) * t));
}

// ----------------------------------------------------------------------------

vec3 CalculateSurfaceNormal(const vec3& p)
{
    const float H = 0.001f;
    const float dx = Density_Func(p + vec3(H, 0.f, 0.f)) - Density_Func(p - vec3(H, 0.f, 0.f));
    const float dy = Density_Func(p + vec3(0.f, H, 0.f)) - Density_Func(p - vec3(0.f, H, 0.f));
    const float dz = Density_Func(p + vec3(0.f, 0.f, H)) - Density_Func(p - vec3(0.f, 0.f, H));

    return glm::normalize(vec3(dx, dy, dz));
}
// ----------------------------------------------------------------------------
RelativePosition triangleRelativePosition(const DefaultMesh &mesh, const DefaultMesh::FaceHandle &faceHandle, glm::vec3 min, float size) {
    //trace("triangle relative position");
    vec3 max = min + size*CHILD_MIN_OFFSETS[7];

    DefaultMesh::FaceVertexIter fv_it = mesh.cfv_iter(faceHandle);
    DefaultMesh::VertexHandle verticesHandle[3];
    DefaultMesh::Point a, b, c;
    int index = 0;
    while (fv_it.is_valid())
    {
        verticesHandle[index++] = *fv_it;
        ++fv_it;
    }
    a = mesh.point(verticesHandle[0]);
    b = mesh.point(verticesHandle[1]);
    c = mesh.point(verticesHandle[2]);

    if (((a[0] > max.x) && (b[0] > max.x) && (c[0] > max.x))
        || ((a[0] < min.x) && (b[0] < min.x) && (c[0] < min.x))){
        return OUTSIDE;
    }
    if (((a[1] > max.y) && (b[1] > max.y) && (c[1] > max.y))
        || ((a[1] < min.y) && (b[1] < min.y) && (c[1] < min.y))){
        return OUTSIDE;
    }
    if (((a[2] > max.z) && (b[2] > max.z) && (c[2] > max.z))
        || ((a[2] < min.z) && (b[2] < min.z) && (c[2] < min.z))){
        return OUTSIDE;
    }

    int numVerticesInside = 0;
    if (vertexRelativePosition(mesh, verticesHandle[0], min, size) == INSIDE){
        ++numVerticesInside;
    }
    if (vertexRelativePosition(mesh, verticesHandle[1], min, size) == INSIDE){
        ++numVerticesInside;
    }
    if (vertexRelativePosition(mesh, verticesHandle[2], min, size) == INSIDE){
        ++numVerticesInside;
    }

    if (numVerticesInside == 3) return INSIDE;

    return CROSSING;
}

// -------------------------------------------------------------------------------
RelativePosition vertexRelativePosition(const DefaultMesh &mesh, const DefaultMesh::VertexHandle &vertexHandle, glm::vec3 min, float size) {
    /*We ignore the case when the vertex is exactly on the cell edge (we consider it outside) */

    vec3 maxvertex = min + (CHILD_MIN_OFFSETS[7] * size);
    //trace("vertex relative position");
    DefaultMesh::Point vertex = mesh.point(vertexHandle);
    //trace("xabu?");
    if ((vertex[0] > min.x) && (vertex[0] < maxvertex.x)
        && (vertex[1] > min.y) && (vertex[1] < maxvertex.y)
        && (vertex[2] > min.z) && (vertex[2] < maxvertex.z)) {
        return INSIDE;
    }
    return OUTSIDE;
}

// -------------------------------------------------------------------------------
RelativePosition halfedgeRelativePosition(const DefaultMesh &mesh, const DefaultMesh::HalfedgeHandle &halfedgeHandle, glm::vec3 min, float size) {
    int numVerticesInside = 0;
    if (vertexRelativePosition(mesh, mesh.to_vertex_handle(halfedgeHandle), min, size) == INSIDE){
        ++numVerticesInside;
    }
    if (vertexRelativePosition(mesh, mesh.from_vertex_handle(halfedgeHandle), min, size) == INSIDE){
        ++numVerticesInside;
    }
    if (numVerticesInside == 2) return INSIDE;
    if (numVerticesInside == 0) return OUTSIDE;
    return CROSSING;
}
// -------------------------------------------------------------------------------
RelativePosition triangleRelativePosition(const Vertex &a, const Vertex &b, const Vertex &c, glm::vec3 min, float size) {
    vec3 max = min + size*CHILD_MIN_OFFSETS[7];
    if (((a.position.x > max.x) && (b.position.x > max.x) && (c.position.x > max.x))
        || ((a.position.x < min.x) && (b.position.x < min.x) && (c.position.x < min.x))){
        return OUTSIDE;
    }
    if (((a.position.y > max.y) && (b.position.y > max.y) && (c.position.y > max.y))
        || ((a.position.y < min.y) && (b.position.y < min.y) && (c.position.y < min.y))){
        return OUTSIDE;
    }
    if (((a.position.z > max.z) && (b.position.z > max.z) && (c.position.z > max.z))
        || ((a.position.z < min.z) && (b.position.z < min.z) && (c.position.z < min.z))){
        return OUTSIDE;
    }

    int numVerticesInside = 0;
    if (vertexRelativePosition(a, min, size) == INSIDE){
        ++numVerticesInside;
    }
    if (vertexRelativePosition(b, min, size) == INSIDE){
        ++numVerticesInside;
    }
    if (vertexRelativePosition(c, min, size) == INSIDE){
        ++numVerticesInside;
    }

    if (numVerticesInside == 3) return INSIDE;


    return CROSSING;
}

// -------------------------------------------------------------------------------

RelativePosition vertexRelativePosition(const Vertex &vertex, glm::vec3 min, float size) {
    /*We ignore the case when the vertex is exactly on the cell edge (we consider it outside) */

    vec3 maxvertex = min + (CHILD_MIN_OFFSETS[7] * size);
    if ((vertex.position.x > min.x) && (vertex.position.x < maxvertex.x)
        && (vertex.position.y > min.y) && (vertex.position.y < maxvertex.y)
        && (vertex.position.z > min.z) && (vertex.position.z < maxvertex.z)) {
        return INSIDE;
    }
    return OUTSIDE;
}

// -------------------------------------------------------------------------------
glm::vec3 openmesh_to_glm(const OpenMesh::VectorT<float, 3> om_vec)
{
    return glm::make_vec3(&om_vec[0]);
}

// -------------------------------------------------------------------------------

bool ray_box_overlap(const glm::vec3 origin, const glm::vec3 dest, const glm::vec3 min, const float size){
    //TODO: implement
    return true;
}

// -------------------------------------------------------------------------------
int compute_nearmost_index(glm::vec3 pivot, std::vector<glm::vec3> intersections)
{
    int index = 0;
    Real mindistance = 99999999.9;
    for (int i = 0; i < intersections.size(); ++i) {
        Real sqrdist = glm::distance2(pivot, intersections[i]);
        if ( sqrdist < mindistance)
        {
            mindistance = sqrdist;
            index = i;
        }
    }
    return index;
}

// -------------------------------------------------------------------------------
void print_point(glm::vec3 point){
    std::cout << "(" << point.x << ", " << point.y << ", " << point.z << ") ";
}

// -------------------------------------------------------------------------------
void print_points(std::vector<glm::vec3> points){
    for (std::vector<glm::vec3>::iterator p_it = points.begin(); p_it != points.end(); ++p_it) {
        print_point(*p_it);
    }
    std::cout << std::endl;
}
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

#include "density.h"
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
void updateSignsArray(int *vecsigns, int size)
{
    bool checksigns = true;
    while(checksigns) {
        checksigns = false;
        for (size_t i = 0; i < size; ++i) {
            if (vecsigns[i] == MATERIAL_UNKNOWN) {
                checksigns = true;
                for (int j = 0; j < 3; ++j) {
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

void write_Ply(std::vector<glm::vec3> &vertices, std::vector<int> &faces, std::string filename)
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
        << "element face " << faces.size()/3 << std::endl
        << "property list uchar int vertex_index" << std::endl
        << "end_header" << std::endl;
        for ( std::vector<glm::vec3>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) {
            outfile << v_it->x << " " << v_it->y << " " << v_it->z << std::endl;
        }

        for ( std::vector<int>::iterator f_it = faces.begin(); f_it != faces.end(); f_it+=3) {
            outfile << "3" << " " << *f_it << " " << *(f_it+1) << " " << *(f_it+2) << std::endl;
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

//// DEBUG ------------------------------------------------------
//if (node->innerFaces.size() > 0)
//{
//std::cout << "On Height Size: " << node->depth << " " << node->size << " Inner Before: " << node->innerFaces.size() << std::endl;
//}
//if (node->crossingFaces.size() > 0)
//{
//std::cout << "On Height: " << node->depth << " " << node->size << " Crossing Before: " << node->crossingFaces.size() << std::endl;
//}
//// DEBUG ------------------------------------------------------ //
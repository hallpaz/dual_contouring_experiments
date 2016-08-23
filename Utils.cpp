#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Utils.h"

using glm::vec3;

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
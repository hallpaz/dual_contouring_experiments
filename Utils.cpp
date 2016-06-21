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

void write_OFF(std::vector<Vertex> &vertices, std::vector<Triangle> &faces, std::string filename)
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

float read_OFF(std::vector<Vertex> &vertices, std::vector<Triangle> &faces, std::string filename) {
    std::ifstream inputfile;
    inputfile.open(filename.c_str(), std::ios::out);
    float mind = 999999.9, maxd = -999999.9;
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
            if (x < mind) mind = x;
            if (y < mind) mind = y;
            if (z < mind) mind = z;
            if (x > maxd) maxd = x;
            if (y > maxd) maxd = y;
            if (z > maxd) maxd = z;
            vertices.push_back(Vertex(glm::vec3(x, y, z)));
        }

        for (int i = 0; i < numFaces; ++i) {
            int n, a, b, c;
            inputfile >> n >> a >> b >> c;
            faces.push_back(Triangle{a, b, c});
        }
        inputfile.close();
    }
    float size = fabs(maxd) > fabs(mind) ? ceilf(fabs(maxd)) : ceilf(fabs(mind));
    std::cout << "Bounding box half size: " << size << std::endl;
    return 2*size;
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
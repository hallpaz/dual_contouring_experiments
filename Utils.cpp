#include <fstream>
#include <iostream>
#include "Utils.h"


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
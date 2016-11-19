//
// Created by hallpaz on 9/1/16.
//

#include "NormalsEstimator.h"
#include <iostream>

using namespace OpenMesh;

void NormalsEstimator::update_model_normals(std::string model_path, std::string output_path, Quality q)
{
    DefaultMesh mesh;
    OpenMesh::IO::Options options;
    if ( ! OpenMesh::IO::read_mesh(mesh, model_path, options) )
    {
        std::cerr << "Error: Cannot read mesh from " << model_path << std::endl;
        return;
    }
    if (!options.check(OpenMesh::IO::Options::FaceNormal))
    {
        mesh.request_face_normals();
        if (!mesh.has_face_normals())
        {
            std::cerr << "Error: Face normals are not available " << std::endl;
            return;
        }
        mesh.update_face_normals();
    }

    mesh.request_vertex_normals();
    if (!mesh.vertex_normals())
    {
        std::cerr << "Error: Vertex normals are not available" << std::endl;
    }

    if (q == COARSE){
        // fast mode
        //mesh.update_vertex_normals();
    }
    else
    {
        DefaultMesh::VertexIter  v_it(mesh.vertices_begin()), v_end(mesh.vertices_end());
        for (; v_it!=v_end; ++v_it)
        {
            DefaultMesh::Normal n;
            mesh.calc_vertex_normal_correct(*v_it, n);
            DefaultMesh::Scalar norm = n.length();
            if (norm != 0.0)
                n *= (DefaultMesh::Scalar(1.0)/norm);
            mesh.set_normal(*v_it, n);
        }
    }
    options += OpenMesh::IO::Options::VertexNormal;
    if ( !OpenMesh::IO::write_mesh(mesh, output_path, options) )
    {
        std::cerr << "Cannot write mesh to file " << output_path << std::endl;
        return;
    }
    else
    {
        std::cout << "Mesh written to " << output_path << std::endl;
    }
}

void NormalsEstimator::compute_better_normals(DefaultMesh &mesh)
{
    mesh.request_face_normals();
    if (!mesh.has_face_normals())
    {
        std::cerr << "Error: Face normals are not available " << std::endl;
        return;
    }
    mesh.update_face_normals();

    mesh.request_vertex_normals();
    if (!mesh.vertex_normals())
    {
        std::cerr << "Error: Vertex normals are not available" << std::endl;
    }


    DefaultMesh::VertexIter v_it(mesh.vertices_begin()), v_end(mesh.vertices_end());
    for (; v_it!=v_end; ++v_it)
    {
        DefaultMesh::Normal n;
        mesh.calc_vertex_normal_correct(*v_it, n);
        DefaultMesh::Scalar norm = n.length();
        if (norm != 0.0)
            n *= (DefaultMesh::Scalar(1.0)/norm);
        mesh.set_normal(*v_it, n);
    }
}
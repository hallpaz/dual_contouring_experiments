cmake_minimum_required(VERSION 3.5)
project(dual_contouring_sandbox)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set(SOURCE_FILES main.cpp old/octree.cpp old/octree.h density.cpp density.h old/qef.cpp old/qef.h old/svd.cpp old/svd.h Utils.cpp Utils.h DataStructures.h old/NormalsEstimator.cpp old/NormalsEstimator.h Contouring.cpp Contouring.h Constants.h old/Tests.cpp old/Tests.h Reconstruction.h Reconstruction.cpp Octree.cpp Octree.h)
add_executable(dual_contouring_sandbox ${SOURCE_FILES})

include(FindOpenMesh.cmake)

include_directories(${OPENMESH_INCLUDE_DIR})

target_link_libraries(dual_contouring_sandbox ${OPENMESH_LIBRARIES})

cmake_minimum_required(VERSION 3.5)
project(dual_contouring_experiments)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set(SOURCE_FILES main.cpp octree.cpp octree.h density.cpp density.h qef.cpp qef.h svd.cpp svd.h Utils.cpp Utils.h DataStructures.h NormalsEstimator.cpp NormalsEstimator.h Contouring.cpp Contouring.h Constants.h Tests.cpp Tests.h)
add_executable(dual_contouring_experiments ${SOURCE_FILES})

include(FindOpenMesh.cmake)

include_directories(${OPENMESH_INCLUDE_DIR})

target_link_libraries(dual_contouring_experiments ${OPENMESH_LIBRARIES})
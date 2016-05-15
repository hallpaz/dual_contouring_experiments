cmake_minimum_required(VERSION 3.5)
project(dual_contouring_experiments)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp octree.cpp octree.h density.cpp density.h qef.cpp qef.h svd.cpp svd.h Utils.cpp Utils.h dataStructures.h)
add_executable(dual_contouring_experiments ${SOURCE_FILES})
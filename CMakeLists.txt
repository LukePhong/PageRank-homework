cmake_minimum_required(VERSION 3.12)
project(pagerank)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp")

add_executable(pagerank main.cpp Pagerank.cpp Pagerank.h)
add_compile_options(-O3)
add_compile_options(-march=corei7-avx)
add_compile_options(-Wall)
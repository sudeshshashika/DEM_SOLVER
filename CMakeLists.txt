cmake_minimum_required(VERSION 3.16)

project(DEM_SOLVER LANGUAGES CXX CUDA)
           
#enable_language(CUDA)
add_executable(DEM_SOLVER DEM.cu)


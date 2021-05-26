cmake_minimum_required(VERSION 3.14)
project(Density)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory("${PROJECT_SOURCE_DIR}/HWC" "${PROJECT_SOURCE_DIR}/HWC/cmake-build-debug")

add_subdirectory("${PROJECT_SOURCE_DIR}/ESM" "${PROJECT_SOURCE_DIR}/ESM/cmake-build-debug")

#add_executable(dense.out densestHyperSubGraph-v2-MultiplyEdges-MultiInstance-newDelete-clearpending.cpp)

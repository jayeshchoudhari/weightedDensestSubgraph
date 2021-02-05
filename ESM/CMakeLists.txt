cmake_minimum_required(VERSION 3.11)
project(Density)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ExecAdd.out ExecAdd.cpp DynGraph.cpp UDynGraph.cpp DynGraphUtils.cpp DSAlgs.cpp DynDSAlg.cpp DynDSAlgAddRem.cpp GraphScheduler.cpp Stats.cpp)
add_executable(ExecAddRem.out ExecAddRem.cpp DynGraph.cpp UDynGraph.cpp DynGraphUtils.cpp DSAlgs.cpp DynDSAlg.cpp DynDSAlgAddRem.cpp GraphScheduler.cpp Stats.cpp)
cmake_minimum_required(VERSION 3.11)
project(Density)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(add.out Add.cpp Hypergraph.cpp ApproximateDS.cpp GraphScheduler.cpp Utility.cpp Stats.cpp ApproximateDSFullyDyn.cpp Decompose.cpp)
add_executable(addRemove.out AddRemove.cpp Hypergraph.cpp ApproximateDS.cpp GraphScheduler.cpp Utility.cpp Stats.cpp ApproximateDSFullyDyn.cpp Decompose.cpp)
add_executable(Approximate.out Approximate.cpp Hypergraph.cpp ApproximateDS.cpp GraphScheduler.cpp Utility.cpp Stats.cpp ApproximateDSFullyDyn.cpp Decompose.cpp)
add_executable(test.out test.cpp)
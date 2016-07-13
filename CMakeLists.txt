cmake_minimum_required(VERSION 3.5)
project(HyperplaneFinderD4)

include_directories(include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES
        src/main.cpp
        src/Line.cpp
        src/CombinationGenerator.cpp
        src/PointsGeometry.cpp)

add_executable(HyperplaneFinderD4 ${SOURCE_FILES})

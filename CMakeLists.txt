cmake_minimum_required(VERSION 3.0)
project(HyperplaneFinderD4)

include_directories(include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -pedantic -pedantic-errors -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wstrict-overflow=5 -Wredundant-decls -Wshadow -Wsign-conversion -Wundef -Wfatal-errors -Wvector-operation-performance -Winvalid-pch -Wconversion -Wdouble-promotion -Wunreachable-code -Winline -Wfloat-equal -Wold-style-cast -Weffc++ -Woverloaded-virtual -Wstrict-null-sentinel -Wsign-promo -Wnoexcept -Wctor-dtor-privacy -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wsuggest-final-types -Wzero-as-null-pointer-constant -Wuseless-cast")

set(SOURCE_FILES
        src/main.cpp
        src/Line.cpp
        src/CombinationGenerator.cpp
        src/PointsGeometry.cpp)

add_executable(HyperplaneFinderD4 ${SOURCE_FILES})

cmake_minimum_required(VERSION 3.7)
project(HyperplaneFinder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pedantic -pedantic-errors -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wstrict-overflow=5 -Wredundant-decls -Wshadow -Wsign-conversion -Wundef -Wfatal-errors -Wvector-operation-performance -Winvalid-pch -Wconversion -Wdouble-promotion -Wunreachable-code -Wfloat-equal -Wold-style-cast -Weffc++ -Woverloaded-virtual -Wstrict-null-sentinel -Wsign-promo -Wnoexcept -Wctor-dtor-privacy -Wdelete-non-virtual-dtor -Wnon-virtual-dtor -Wsuggest-final-types -Wzero-as-null-pointer-constant -Wuseless-cast")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
include_directories(src)

set(SOURCE_FILES src/main.cpp)
add_executable(HyperplaneFinder ${SOURCE_FILES})

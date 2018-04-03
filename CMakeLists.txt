cmake_minimum_required(VERSION 3.10) # For CXX_STANDARD 17 property on Visual Studio

project(HyperplaneFinder)
enable_language(CXX)

include(cmake/common.cmake)

set(CMAKE_DEBUG_POSTFIX _d)

# By default build in Release mode
if( NOT CMAKE_BUILD_TYPE )
	set(CMAKE_BUILD_TYPE "Release")
endif()

get_filename_component(root_dir   ${CMAKE_SOURCE_DIR}   REALPATH)
get_filename_component(src_dir    ${root_dir}/src       REALPATH)

get_files(
  files "${src_dir}"
  OPTIONS recurse
)

make_target(
  HyperplaneFinder "HyperplaneFinder" ${files}
  INCLUDES "${src_dir}"
  OPTIONS cxx executable
)
set_property(TARGET HyperplaneFinder PROPERTY CXX_STANDARD 17)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	target_add_flag(HyperplaneFinder "-march=native" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-msse4.2" DEBUG RELEASE RELWITHDEBINFO)
endif()

# Debug
#set(CMAKE_VERBOSE_MAKEFILE 1)

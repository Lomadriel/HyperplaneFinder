cmake_minimum_required(VERSION 3.10) # For CXX_STANDARD 17 property on Visual Studio

project(HyperplaneFinder)
enable_language(CXX)

include(cmake/common.cmake)
include(cmake/json.cmake)
include(cmake/inja.cmake)

set(CMAKE_DEBUG_POSTFIX _d)

# By default build in Release mode
if( NOT CMAKE_BUILD_TYPE )
	set(CMAKE_BUILD_TYPE "Release")
endif()

get_filename_component(root_dir        ${CMAKE_SOURCE_DIR}     ABSOLUTE)
get_filename_component(src_dir         ${root_dir}/src         ABSOLUTE)
get_filename_component(templates_dir   ${root_dir}/templates   ABSOLUTE)

get_files(
  files "${src_dir}"
  OPTIONS recurse
)

make_target(
  HyperplaneFinder "HyperplaneFinder" ${files}
  INCLUDES "${src_dir}"
  EXT_INCLUDES "${JSON_INCLUDE_DIR}" "${INJA_INCLUDE_DIR}"
  OPTIONS cxx executable
)
if(NOT MSVC)
	target_link_libraries(HyperplaneFinder stdc++fs)
endif()
set_property(TARGET HyperplaneFinder PROPERTY CXX_STANDARD 17)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	# https://gcc.gnu.org/onlinedocs/gcc-4.5.3/gcc/i386-and-x86_002d64-Options.html
	target_add_flag(HyperplaneFinder "-march=native" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-mtune=native" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-msse4.2" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-mavx" DEBUG RELEASE RELWITHDEBINFO)

	# https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html#Optimize-Options
	target_add_flag(HyperplaneFinder "-fmodulo-sched" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fmodulo-sched-allow-regmoves" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fgcse-sm" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fgcse-las" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fgcse-after-reload" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fdevirtualize-speculatively" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fdevirtualize-at-ltrans" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-floop-nest-optimize" DEBUG RELEASE RELWITHDEBINFO) # experimental
	target_add_flag(HyperplaneFinder "-floop-parallelize-all" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-fsplit-loops" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-funswitch-loops" DEBUG RELEASE RELWITHDEBINFO)
	target_add_flag(HyperplaneFinder "-floop-unroll-and-jam" DEBUG RELEASE RELWITHDEBINFO)
endif()

configure_folder(${templates_dir}   "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/templates")

# Debug
#set(CMAKE_VERBOSE_MAKEFILE 1)

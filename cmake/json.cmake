message(STATUS "Configuring JSON")

# Config
set(JSON_NAME      "json")
set(JSON_VERSION   "3.1.2")

get_filename_component(JSON_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/extlibs/${JSON_NAME}-${JSON_VERSION}/include ABSOLUTE)

# Message
message("> include: ${JSON_INCLUDE_DIR}")

message(STATUS "Configuring INJA")

# Config
set(INJA_NAME      "inja")
set(INJA_VERSION   "1.0.0")

get_filename_component(INJA_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/extlibs/${INJA_NAME}-${INJA_VERSION}/include REALPATH)

# Message
message("> include: ${INJA_INCLUDE_DIR}")

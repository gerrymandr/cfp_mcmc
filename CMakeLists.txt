cmake_minimum_required (VERSION 2.8)
project (chain)

#version number
set (chain_VERSION_MAJOR 1)
set (chain_VERSION_MINOR 0)

set (CMAKE_CXX_STANDARD 11)

#indicate the entry point for the executable
FILE(GLOB SRC_FILES "${CMAKE_SOURCE_DIR}/src/*")
add_executable (chain ${SRC_FILES})

install (TARGETS chain DESTINATION bin)

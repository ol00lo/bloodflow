cmake_minimum_required(VERSION 3.8)

Project(bflow)

include_directories("${CMAKE_SOURCE_DIR}/extern/catch2")
include_directories("${CMAKE_SOURCE_DIR}/src/bflow")

add_subdirectory(src)


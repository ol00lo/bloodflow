cmake_minimum_required(VERSION 3.0)

Project(bflow)

include_directories("${CMAKE_SOURCE_DIR}/extern/catch2")

add_subdirectory(src)

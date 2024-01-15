cmake_minimum_required(VERSION 3.8)

Project(bflow)

include_directories("${CMAKE_SOURCE_DIR}/extern/catch2")

add_subdirectory(src)

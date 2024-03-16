cmake_minimum_required(VERSION 3.8)

Project(bflow)

include_directories("${CMAKE_SOURCE_DIR}/extern/catch2")
include_directories("${CMAKE_SOURCE_DIR}/extern")
include_directories("${CMAKE_SOURCE_DIR}/src")

set(BOOST_ROOT ${CMAKE_SOURCE_DIR}/../../)
find_package(Boost REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(src)


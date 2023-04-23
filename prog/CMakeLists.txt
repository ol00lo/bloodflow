cmake_minimum_required(VERSION 3.0)

Project(bloodflow)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

# build all binaries to build/bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ------ amgcl
FetchContent_Declare(
	amgcl
	GIT_REPOSITORY https://github.com/ddemidov/amgcl
	GIT_TAG 1.4.2
)
FetchContent_GetProperties(amgcl)
if(NOT amgcl_POPULATED)
	FetchContent_Populate(amgcl)
endif()
include_directories(${amgcl_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/../../boost_1_80_0)
set(BLOODFLOWLIB_TARGET bloodflow)

add_subdirectory(lib)
add_subdirectory(tests)

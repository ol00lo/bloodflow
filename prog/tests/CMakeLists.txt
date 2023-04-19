# Fetch catch2 unit testing framework
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.5
)

FetchContent_GetProperties(catch)
if(NOT catch_POPULATED)
  FetchContent_Populate(catch)
  add_subdirectory(${catch_SOURCE_DIR} ${catch_BINARY_DIR})
endif()

set(CATCH2_FRAMEWORK "Catch2::Catch2")

set(SRC
	transport_test.cpp
	bloodflow_test.cpp
	model_test.cpp
)
set(HEADERS
	testutils.hpp)

add_executable(bloodflow_test ${SRC})
target_link_libraries(bloodflow_test PRIVATE ${BLOODFLOWLIB_TARGET} ${CATCH2_FRAMEWORK})

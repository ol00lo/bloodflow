
set (SRC
	tests/bflow_test.cpp
	tests/graph_test.cpp
	tests/grid_test.cpp
	tests/vtk_test.cpp
	tests/test_utils.cpp
	bflow/vessel_graph.cpp
	bflow/graph_grid.cpp
	bflow/vtk.cpp
	bflow/matrix.cpp
	tests/matrix_test.cpp
	bflow/debug/printer.cpp
)
set(HEADERS
	bflow/vessel_graph.hpp
	bflow/graph_grid.hpp
	bflow/vtk.hpp
	tests/test_utils.hpp
	bflow/matrix.hpp
	bflow/debug/printer.hpp
)
add_executable(bflow_test ${SRC} ${HEADERS})

set (SRC
	tests/bflow_test.cpp
	tests/graph_test.cpp
	tests/grid_test.cpp
	tests/vtk_test.cpp
	bflow/vessel_graph.cpp
	bflow/graph_grid.cpp
	bflow/vtk.cpp
)
set(HEADERS
	bflow/vessel_graph.hpp
	bflow/graph_grid.hpp
	bflow/vtk.hpp
)
add_executable(bflow_test ${SRC} ${HEADERS})
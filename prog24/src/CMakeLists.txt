
set (SRC
	tests/bflow_test.cpp
	tests/graph_test.cpp
	tests/grid_test.cpp
	bflow/vessel_graph.cpp
	bflow/graph_grid.cpp
)
set(HEADERS
	bflow/vessel_graph.hpp
	bflow/graph_grid.hpp
)
add_executable(bflow_test ${SRC} ${HEADERS})
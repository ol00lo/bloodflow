
set (SRC
	tests/bflow_test.cpp
	tests/graph_test.cpp
	bflow/vessel_graph.cpp
)
set(HEADERS
	bflow/vessel_graph.hpp
)
add_executable(bflow_test ${SRC} ${HEADERS})
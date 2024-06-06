
set (SRC
	tests/bflow_test.cpp
	tests/graph_test.cpp
	tests/grid_test.cpp
	tests/vtk_test.cpp
	tests/test_utils.cpp
	tests/matrix_test.cpp
	tests/matrix_solver_test.cpp
	tests/dg_test.cpp
	tests/vessels_dg_test.cpp
	tests/nonlinsolver.cpp
	bflow/vessel_graph.cpp
	bflow/graph_grid.cpp
	bflow/vtk.cpp
	bflow/matrix.cpp
	bflow/debug/printer.cpp
	bflow/matrix_solver.cpp
	bflow/time_series_writer.cpp
)
set(HEADERS
	tests/test_utils.hpp
	tests/nonlinsolver.hpp
	bflow/vessel_graph.hpp
	bflow/graph_grid.hpp
	bflow/vtk.hpp
	bflow/matrix.hpp
	bflow/debug/printer.hpp
	bflow/matrix_solver.hpp
	bflow/time_series_writer.hpp
	bflow/macros.hpp
)
add_executable(bflow_test ${SRC} ${HEADERS})
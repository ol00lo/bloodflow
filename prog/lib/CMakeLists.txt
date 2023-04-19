SET (HEADERS
	common.hpp
	grid.hpp
	debug.hpp
	single_vessel.hpp
	slae/csrmat.hpp
	slae/matrix_solver.hpp
	assem/transport.hpp
	assem/misc.hpp
	logger/imonitor.hpp
	logger/logger.hpp
	logger/monitor_trigger.hpp
	logger/vtk_monitor.hpp
)

SET (SRC
	grid.cpp
	debug.cpp
	single_vessel.cpp
	slae/matrix_solver.cpp
	slae/csrmat.cpp
	assem/transport.cpp
	assem/misc.cpp
	logger/logger.cpp
	logger/vtk_monitor.cpp
)


add_library(${BLOODFLOWLIB_TARGET} STATIC ${HEADERS} ${SRC})

#include "common.hpp"
#include "testutils.hpp"
#include "single_vessel.hpp"
#include "grid.hpp"
#include "logger/logger.hpp"
#include "logger/vtk_monitor.hpp"
#include "logger/monitor_trigger.hpp"

using namespace bf;

TEST_CASE("Single vessel, no friction", "[single_vessel][no_friction]"){
	// parameters
	double L = 10; // m
	double A0 = M_PI*1e-4; // m2
	double h = 1.5e-3; // m
	double rho = 1050; // kg/m3
	double mu = 0*1e-3; // Pa*s
	double zeta = 9;
	double E = 4e5; // Pa

	// nondimensional parameters
	double x0 = 1; // m
	double u0 = 1e-6/A0; // m/s
	double q0 = u0*A0;
	double t0 = x0/u0;
	double delta_p = rho*u0*u0;
	double Mp = 4.0/3.0*sqrt(M_PI)*E*h/sqrt(A0)/delta_p;
	double Mf = 2*(zeta+2)*mu*M_PI*x0/(rho*u0*A0);
	std::cout << Mp << " " << Mf << std::endl;

	Grid grid(L/x0, 200);
	auto input_q = [t0, q0](double t){ return 1e-6*exp(-1e4*(t*t0-0.05)*(t*t0-0.05))/q0; };
	SingleVessel problem(grid, Mp, Mf, input_q);

	std::shared_ptr<IMonitor> saver(new VtkMonitor(grid, "single_vessel", {"velocity", "pressure", "area"}));
	std::shared_ptr<IMonitorTrigger> trigger(new MonitorTrigger_TimePeriod(0.05));
	Logger logger;
	logger.add_monitor(saver, trigger);
	
	size_t iter = 0;
	double tau = 1e-6;
	while (problem.time() <= 1000*tau){
		problem.step(tau);
		iter += 1;
		logger.step(iter, problem.time(), problem.solution());
	}
}

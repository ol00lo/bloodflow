
#include "testutils.hpp"
#include <fstream>
#include "assem/misc.hpp"
#include "grid.hpp"
#include "assem/transport.hpp"
#include "debug.hpp"
#include "logger/logger.hpp"
#include "logger/vtk_monitor.hpp"
#include "logger/monitor_trigger.hpp"

using namespace bf;

TEST_CASE("Transport equation, upwind", "[upwind]"){
	auto u0 = [](double x) -> double{
		return (x < -0.2 || x > 0.2) ? 0 : 1;
	};

	auto ufun = [&](double x, double t) -> double{
		return u0(x-t);
	};
	double len = 5.0;
	Grid grid(len, 30*len);

	std::vector<double> u(grid.n_points());
	for (size_t i=0; i<grid.n_points(); ++i){
		u[i] = u0(grid.coo(i));
	}

	std::shared_ptr<ITransport> tran;
	int num_section;
	// SECTION("Upwind"){
	// 	num_section = 1;
	// 	tran.reset(new UpwindTransport());
	// }
	// SECTION("TvdSuperbee"){
	// 	num_section = 2;
	// 	tran.reset(new TvdSuperbeeTransport());
	// }
	// SECTION("TvdMc"){
	// 	num_section = 3;
	// 	tran.reset(new TvdMcTransport());
	// }
	SECTION("VanLeer"){
		num_section = 4;
		tran.reset(new TvdVanLeerTransport());
	}
	// SECTION("MinMod"){
	// 	num_section = 5;
	// 	tran.reset(new TvdMinModTransport());
	// }
	// SECTION("Umist"){
	// 	num_section = 6;
	// 	tran.reset(new TvdUmistTransport());
	// }
	//  SECTION("Ospre"){
	//  	num_section = 7;
	//  	tran.reset(new TvdOspreTransport());
	//  }
	// SECTION("Charm"){
	// 	num_section = 8;
	// 	tran.reset(new TvdCharmTransport());
	// }
	std::vector<double> vel(grid.n_points(), 1);
	tran->set_velocity(vel);

	std::shared_ptr<IMonitor> saver(new VtkMonitor(grid, "u_vanleer"));
	std::shared_ptr<IMonitorTrigger> trigger(new MonitorTrigger_TimePeriod(0.05));
	Logger logger;
	logger.add_monitor(saver, trigger);
	
	double time = 0;
	size_t iter = 0;
	double tau = 0.007;
	double h = grid.len(0);
	double nrm;
	std::ofstream norm_out("norms_vanleer.txt");
	while (time <= 1){
		//fluxes
		std::vector<double> fluxes = tran->compute_fluxes(u);
		double flast = u.back()*vel.back();
		fluxes.push_back(flast);
		time += tau;
		// left boundary conditions
		u[0] = ufun(0, time);
		// next layer
		for (size_t i=1; i<grid.n_points(); ++i){
			size_t cell_left = i-1;
			size_t cell_right = i;
			u[i] -= tau/h * (fluxes[cell_right] - fluxes[cell_left]);
		}


		// compute norm
		nrm = 0;
		for (size_t i=0; i<grid.n_points(); ++i){
			double diff = u[i] - ufun(grid.coo(i), time);
			nrm += (diff*diff);
		}
		nrm = sqrt(nrm/grid.n_points());
		iter += 1;
		norm_out << iter << "\t" << time << "\t" << nrm << std::endl;
		logger.step(iter, time, u);
	}

	switch (num_section){
		case 1: CHECK(nrm == Approx(0.1871625127)); break;
		case 2: CHECK(nrm == Approx(0.0947743)); break;
		case 3: CHECK(nrm == Approx(0.0991519263)); break; 
		case 4: CHECK(nrm == Approx(0.1009686601)); break;
		case 5: CHECK(nrm == Approx(0.1187898205)); break;
		case 6: CHECK(nrm == Approx(0.0999570685)); break;
		case 7: CHECK(nrm == Approx(0.1028867349)); break;
		case 8: CHECK(nrm); break;
		default: _THROW_UNREACHABLE_;
	}
}

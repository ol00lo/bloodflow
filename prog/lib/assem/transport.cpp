#include "transport.hpp"

using namespace bf;

namespace{
double superbee(double r){
	if (r<0)   return 0;
	if (r<0.5) return 2*r;
	if (r<1)   return 1;
	if (r<2)   return r;
	return 2;
}
}

void UpwindTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> UpwindTransport::compute_fluxes(const std::vector<double>& u) const{
	size_t n_cells = u.size() - 1;
	std::vector<double> fluxes(n_cells, 0);
	for (size_t i=0; i<n_cells; ++i){
		size_t left_point = i;
		size_t right_point = i+1;
		double check_v = 0.5*(_velocity[left_point] + _velocity[right_point]);
		if (check_v >= 0){
			fluxes[i] = _velocity[left_point] * u[left_point];
		} else {
			fluxes[i] = _velocity[right_point] * u[right_point];
		}
	}
	return fluxes;
}

void TvdSuperbeeTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdSuperbeeTransport::compute_fluxes(const std::vector<double>& u) const{
	size_t n_points = u.size();
	size_t n_cells = u.size() - 1;

	std::vector<double> r(n_points, 1);
	for (size_t i=1; i<n_points-1; ++i){
		double denum = u[i+1] - u[i];
		if (denum == 0) denum = 1e-8;
		r[i] = (u[i] - u[i-1])/denum;
	}

	std::vector<double> fluxes(n_cells, 0);
	for (size_t i=0; i<n_cells; ++i){
		size_t left_point = i;
		size_t right_point = i+1;
		double flux_high, flux_low, phi;
		double check_v = 0.5*(_velocity[left_point] + _velocity[right_point]);

		flux_high = 0.5*(_velocity[left_point] * u[left_point] + _velocity[right_point] * u[right_point]);
		if (check_v >= 0){
			flux_low = _velocity[left_point] * u[left_point];
			phi = superbee(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = superbee(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}

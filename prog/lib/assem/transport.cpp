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
double mc(double r){
	if (r<=0)   return 0;
	if (r<=1.0/3.0) return 2*r;
	if (r<=3.0)   return (1+r)/2;
	return 2;
}
double vanleer(double r){
	if(r>0) return 2*r/(1+r);
	if(r<=0) return 0;
}
double minmod(double r){
	if (r<=0) return 0;
	if (r<=1) return r;
	return 1;
}
double umist(double r){
	if (r<=0) return 0;
	if (r<=0.2) return 2*r;
	if (r<=3.0/7.0) return 0.25+0.75*r;
	if (r<=7.0/3.0) return 0.75+0.25*r;
	return 2;
}
double ospre(double r){
	if(r>0) return 1.5*(r*r+r)/(r*r+1+r);
	if(r<=0) return 0;
}
double charm(double r){
	if(r<=0) return 0;
	if(r>0) return r*(3*r+1)/(r+1)/(r+1);
	return 3;
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

void TvdMcTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdMcTransport::compute_fluxes(const std::vector<double>& u) const{
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
			phi = mc(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = mc(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}

void TvdVanLeerTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdVanLeerTransport::compute_fluxes(const std::vector<double>& u) const{
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
			phi = vanleer(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = vanleer(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}

void TvdMinModTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdMinModTransport::compute_fluxes(const std::vector<double>& u) const{
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
			phi = minmod(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = minmod(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}

void TvdUmistTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdUmistTransport::compute_fluxes(const std::vector<double>& u) const{
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
			phi = umist(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = umist(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}

void TvdOspreTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdOspreTransport::compute_fluxes(const std::vector<double>& u) const{
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
			phi = ospre(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = ospre(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}

void TvdCharmTransport::set_velocity(const std::vector<double>& vel){
	_velocity = vel;
}

std::vector<double> TvdCharmTransport::compute_fluxes(const std::vector<double>& u) const{
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
			phi = charm(r[left_point]);
		} else {
			flux_low = _velocity[right_point] * u[right_point];
			phi = charm(r[right_point]);
		}
		fluxes[i] = flux_low + phi*(flux_high - flux_low);
	}
	return fluxes;
}
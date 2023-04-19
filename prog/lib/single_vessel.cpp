#include "single_vessel.hpp"
#include "assem/transport.hpp"

using namespace bf;

SingleVessel::SingleVessel(const Grid& grid, double Mp, double Mf, std::function<double(double)> input_q):
	_grid(&grid), _Mp(Mp), _Mf(Mf), _input_q(input_q),
	_pressure(_grid->n_points(), 0),
	_velocity(_grid->n_points(), 0),
	_area(_grid->n_points(), 1){}

void SingleVessel::step(double tau){
	_time += tau;
	double h = _grid->len(0);
	TvdSuperbeeTransport tran;
	tran.set_velocity(_velocity);

	// 1. ==== solve area problem
	//fluxes
	std::vector<double> area_fluxes = tran.compute_fluxes(_area);
	double area_flast = _area.back()*_velocity.back();
	area_fluxes.push_back(area_flast);

	// left boundary condition
	_area[0] = 1;
	// next layer
	for (size_t i=1; i<_grid->n_points(); ++i){
		size_t cell_left = i-1;
		size_t cell_right = i;
		_area[i] -= tau/h * (area_fluxes[cell_right] - area_fluxes[cell_left]);

		// asserts
		if (_area[i] <= 0 || std::isnan(_area[i])){
			_THROW_UNREACHABLE_;
		}
	}

	// 2. === solve pressure problem
	for (size_t i=0; i<_grid->n_points(); ++i){
		_pressure[i] = _Mp*(sqrt(_area[i]) - 1);

		// asserts
		if (std::isnan(_pressure[i])){
			_THROW_UNREACHABLE_;
		}
	}

	// 3. === pressure gradient
	std::vector<double> dpdx(_grid->n_points());
	dpdx[0] = (_pressure[1] - _pressure[0])/h;
	dpdx.back() = (_pressure.back() - _pressure[_grid->n_points()-2])/h;
	for (size_t i=1; i<_grid->n_points()-1; ++i){
		dpdx[i] = (_pressure[i+1] - _pressure[i-1])/(2*h);
	}

	// 4. === solve velocity problem
	//fluxes
	std::vector<double> vel_fluxes = tran.compute_fluxes(_velocity);
	double vel_flast = _velocity.back()*_velocity.back();
	vel_fluxes.push_back(vel_flast);
	// left boundary condition
	_velocity[0] = _input_q(_time);
	std::cout << _input_q(_time) << std::endl;
	// next layer
	for (size_t i=1; i<_grid->n_points(); ++i){
		size_t cell_left = i-1;
		size_t cell_right = i;
		double convection = tau/h * (vel_fluxes[cell_right] - vel_fluxes[cell_left]);
		double coef = 1 + tau*_Mf/_area[i];
		_velocity[i] = (1.0/coef) * (_velocity[i] - convection - tau*dpdx[i]);

		// asserts
		if (std::isnan(_velocity[i])){
			_THROW_UNREACHABLE_;
		}
	}
}

const std::vector<double>& SingleVessel::pressure() const{
	return _pressure;
}

const std::vector<double>& SingleVessel::area() const{
	return _area;
}

const std::vector<double>& SingleVessel::velocity() const{
	return _velocity;
}

std::vector<const std::vector<double>*> SingleVessel::solution() const{
	return {&_velocity, &_pressure, &_area};
}

double SingleVessel::time() const{
	return _time;
}

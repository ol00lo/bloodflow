#ifndef SINGLE_VESSEL_HPP
#define SINGLE_VESSEL_HPP

#include <functional>
#include "grid.hpp"

namespace bf{

class SingleVessel{
public:
	SingleVessel(const Grid& grid, double Mp, double Mf, std::function<double(double)> input_q);
	void step(double tau);

	const std::vector<double>& pressure() const;
	const std::vector<double>& area() const;
	const std::vector<double>& velocity() const;
	std::vector<const std::vector<double>*> solution() const;
	double time() const;
private:
	const Grid* _grid;
	double _Mp;
	double _Mf;
	std::function<double(double)> _input_q;

	// solution data
	std::vector<double> _pressure;
	std::vector<double> _velocity;
	std::vector<double> _area;
	double _time;
};

};

#endif

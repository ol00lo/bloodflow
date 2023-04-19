#ifndef BF_GRID_HPP
#define BF_GRID_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"

namespace bf{

class Grid{
public:
	Grid(double L, double n);
	size_t n_points() const;
	size_t n_cells() const;
	double coo(size_t i) const;
	double len(size_t icell) const;
	const CsrStencil& stencil() const;
	void save_data(const std::vector<double>& v, std::string fn) const;
	void save_data(std::vector<std::string> data_names, std::vector<const std::vector<double>*> data, std::string fn) const;
private:
	const double _L;
	const size_t _n;
	const double _h;
	CsrStencil _stencil;
};

}

#endif

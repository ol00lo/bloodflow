#include "misc.hpp"

using namespace bf;

std::vector<double> bf::assemble_unit(const Grid& grid){
	std::vector<double> ret(grid.n_points(), 0);
	for (size_t i=0; i<grid.n_cells(); ++i){
		double ch = grid.len(i);
		ret[i] += ch/2;
		ret[i+1] += ch/2;
	}

	return ret;
}

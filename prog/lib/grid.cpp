#include "grid.hpp"
#include <fstream>

using namespace bf;

Grid::Grid(double L, double n):_L(L), _n(n), _h(L/n){}

size_t Grid::n_points() const{
	return _n + 1;
}

double Grid::coo(size_t i) const{
	return i*_h;
}

size_t Grid::n_cells() const{
	return _n;
}

double Grid::len(size_t icell) const{
	return _h;
}

const CsrStencil& Grid::stencil() const{
	if (_stencil.n_nonzero()==0){
		CsrStencil* s = const_cast<CsrStencil*>(&_stencil);
		std::vector<size_t> addr {0, 2};
		std::vector<size_t> cols {0, 1};

		for (size_t i=1; i<n_points()-1; ++i){
			cols.push_back(i-1);
			cols.push_back(i);
			cols.push_back(i+1);
			addr.push_back(addr.back()+3);
		}

		cols.push_back(n_points()-2);
		cols.push_back(n_points()-1);
		addr.push_back(addr.back()+2);

		s->set_data(std::move(addr), std::move(cols));
	}

	return _stencil;
}

void Grid::save_data(const std::vector<double>& vec, std::string fn) const{
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Grid1"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET RECTILINEAR_GRID"<<std::endl;
	fs<<"DIMENSIONS " << n_points() << " 1 1" << std::endl;
	fs<<"X_COORDINATES " << n_points() << " double" << std::endl;
	for (size_t i=0; i<n_points(); ++i){
		fs << coo(i) << std::endl;
	}
	fs << "Y_COORDINATES 1 double 0" << std::endl;
	fs << "Z_COORDINATES 1 double 0" << std::endl;
	fs << "POINT_DATA " << n_points() << std::endl;
	fs << "SCALARS data double 1" << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	for (size_t i=0; i<n_points(); ++i){
		fs << vec[i] << std::endl;
	}
	fs.close();
}

void Grid::save_data(std::vector<std::string> data_names, std::vector<const std::vector<double>*> data, std::string fn) const{
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Grid1"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	fs<<"DATASET RECTILINEAR_GRID"<<std::endl;
	fs<<"DIMENSIONS " << n_points() << " 1 1" << std::endl;
	fs<<"X_COORDINATES " << n_points() << " double" << std::endl;
	for (size_t i=0; i<n_points(); ++i){
		fs << coo(i) << std::endl;
	}
	fs << "Y_COORDINATES 1 double 0" << std::endl;
	fs << "Z_COORDINATES 1 double 0" << std::endl;
	fs << "POINT_DATA " << n_points() << std::endl;
	for (size_t idata=0; idata<data_names.size(); ++idata){
		fs << "SCALARS " << data_names[idata] << " double 1" << std::endl;
		fs << "LOOKUP_TABLE default" << std::endl;
		for (size_t i=0; i<n_points(); ++i){
			fs << (*data[idata])[i] << std::endl;
		}
	}
	fs.close();
}

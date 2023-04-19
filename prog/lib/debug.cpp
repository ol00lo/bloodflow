#include "debug.hpp"

using namespace bf;

#include <iomanip>
#include <iostream>
#include "debug.hpp"


void dbg::print(const CsrStencil& mat, const std::vector<double>& vals){
	const int ndigits = 8;

	auto& cols = mat.cols();
	auto& addr = mat.addr();
	int ncols = *std::max_element(cols.begin(), cols.end()) + 1;

	std::cout << std::setw(3) << " " << "  ";
	for (int icol=0; icol<ncols; ++icol){
		std::cout << std::setw(ndigits) << icol << " ";
	}
	std::cout << std::endl;
	std::cout << std::setw(3) << " " << "  ";
	std::string cap(ndigits, '-');
	for (int icol=0; icol<ncols; ++icol){
		std::cout << cap << " ";
	}
	std::cout << std::endl;

	for (size_t irow=0; irow<mat.n_rows(); ++irow){
		std::vector<double> row_vals(ncols, std::numeric_limits<double>::max());
		for (size_t a=addr[irow]; a<addr[irow+1]; ++a){
			row_vals[cols[a]] = vals[a];
		}

		std::cout << std::setw(3) << irow << "| ";
		for (double v: row_vals){
			if (v == std::numeric_limits<double>::max()){
				std::cout << std::setw(ndigits) << "*";
			} else {
				std::cout << std::setw(ndigits) << v;
			}
			std::cout << " ";
		}
		std::cout << std::endl;
	}
}

void dbg::print(const CsrMatrix& mat){
	dbg::print(mat, mat.vals());
}

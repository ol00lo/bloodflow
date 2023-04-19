#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"

namespace bf{

// solve SLAE using amgc method
class AmgcMatrixSolver{
public:
	AmgcMatrixSolver(size_t maxit=1000, double eps=1e-8);
	~AmgcMatrixSolver();

	void set_matrix(const CsrMatrix& mat);
	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values);

	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const;
private:
	size_t _maxit;
	double _tolerance;

	struct Impl;
	std::unique_ptr<Impl> _pimpl;
};

}

#endif

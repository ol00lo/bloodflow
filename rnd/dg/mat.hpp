#ifndef __DG_MAT_CSRMAT_HPP__
#define __DG_MAT_CSRMAT_HPP__

#include "common.hpp"

class CsrStencil{
public:
	virtual ~CsrStencil() = default;
	void set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols);
	void set_stencil(const std::vector<std::set<size_t>>& stencil_set);

	const std::vector<size_t>& addr() const;
	const std::vector<size_t>& cols() const;
	virtual void validate() const;

	size_t get_address(size_t irow, size_t icol) const;

	size_t n_rows() const;
	size_t n_nonzeros() const;
	bool is_in_stencil(size_t irow, size_t icol) const;
private:
	std::vector<size_t> _addr = {0};
	std::vector<size_t> _cols;
};

class CsrMatrix: public CsrStencil{
public:
	CsrMatrix() = default;
	CsrMatrix(const CsrStencil& stencil): CsrStencil(stencil), _vals(stencil.n_nonzeros(), 0.0){}

	void set_values(std::vector<double>&& vals);
	const std::vector<double>& vals() const;
	std::vector<double>& vals();
	void set_unit_row(size_t irow);
	void validate() const;
	double value(size_t irow, size_t icol) const;
	std::vector<double> mult_vec(const std::vector<double>& u) const;
	double mult_vec(size_t irow, const std::vector<double>& u) const;
private:
	std::vector<double> _vals;
};


class AmgcMatrixSolver{
public:
	AmgcMatrixSolver(int maxit=1000, double eps=1e-8);
	AmgcMatrixSolver(std::initializer_list<std::pair<std::string, std::string>> amgc_params);
	~AmgcMatrixSolver();
	void set_matrix(const CsrMatrix& mat);
	void set_matrix(const CsrStencil& mat_stencil, const std::vector<double>& mat_values);

	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const;

	static void solve_slae(const CsrMatrix& mat, const std::vector<double>& rhs, std::vector<double>& x, int maxit=1000, double eps=1e-8);
private:
	int _maxit;
	double _tolerance;
	std::map<std::string, std::string> _params;

	struct Impl;
	std::unique_ptr<Impl> _pimpl;
};

std::ostream& operator<<(std::ostream& os, const CsrMatrix& mat);

#endif

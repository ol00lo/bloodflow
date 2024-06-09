#include "mat.hpp"

#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>

void CsrStencil::set_stencil(std::vector<size_t>&& addr, std::vector<size_t>&& cols){
	_addr = std::move(addr);
	_cols = std::move(cols);
}

void CsrStencil::set_stencil(const std::vector<std::set<size_t>>& stencil_set){
	_addr = std::vector<size_t>(1, 0);
	_cols.clear();

	for (size_t irow=0; irow<stencil_set.size(); ++irow){
		const std::set<size_t>& cols = stencil_set[irow];
		_addr.push_back(_addr.back() + cols.size());
		for (size_t col: cols){
			_cols.push_back(col);
		}
	}
}

size_t CsrStencil::n_nonzeros() const{
	return _cols.size();
}

size_t CsrStencil::n_rows() const{
	return _addr.size()-1;
}

const std::vector<size_t>& CsrStencil::addr() const{
	return _addr;
}

const std::vector<size_t>& CsrStencil::cols() const{
	return _cols;
}

void CsrStencil::validate() const{
	// sizes
	if (_addr.size() < 1){
		throw std::runtime_error("addr array should have more then zero entries");
	}
	if (_cols.size() != _addr.back()){
		throw std::runtime_error("cols array size should match last addr entry");
	}
	// non-decreasing
	for (size_t i=1; i<_addr.size(); ++i){
		if (_addr[i] < _addr[i-1]){
			throw std::runtime_error("addr array should be non-decreasing");
		}
	}
}

bool CsrStencil::is_in_stencil(size_t irow, size_t icol) const{
	size_t start = _addr.at(irow);
	size_t end = _addr.at(irow+1);
	for (size_t i=start; i<end; ++i){
		if (_cols.at(i) == icol){
			return true;
		}
	}
	return false;
}

size_t CsrStencil::get_address(size_t irow, size_t icol) const{
	std::vector<size_t>::const_iterator it_start = _cols.begin() + _addr.at(irow);
	std::vector<size_t>::const_iterator it_end = _cols.begin() + _addr.at(irow+1);
	auto fnd = std::lower_bound(it_start, it_end, icol);
	if (fnd != it_end && *fnd == icol){
		size_t a = fnd - _cols.begin();
		return a;
	}
	return INVALID_INDEX;
}

void CsrMatrix::set_values(std::vector<double>&& vals){
	_vals = std::move(vals);
}

const std::vector<double>& CsrMatrix::vals() const{
	return _vals;
}

std::vector<double>& CsrMatrix::vals(){
	return _vals;
}

void CsrMatrix::validate() const{
	CsrStencil::validate();

	// values size
	if (_vals.size() != n_nonzeros()){
		throw std::runtime_error("values array should have same size as the columns arrays");
	}
}

double CsrMatrix::value(size_t irow, size_t icol) const{
	size_t a = get_address(irow, icol);
	if (a != INVALID_INDEX){
		return _vals[a];
	} else {
		return 0.0;
	}
}

std::vector<double> CsrMatrix::mult_vec(const std::vector<double>& u) const{
	const std::vector<size_t>& a = addr();
	const std::vector<size_t>& c = cols();
	const std::vector<double>& v = vals();

	std::vector<double> ret(n_rows(), 0);
	for (size_t irow=0; irow<n_rows(); ++irow){
		size_t start = a[irow];
		size_t end = a[irow+1];
		for (size_t i=start; i<end; ++i){
			ret[irow] += v[i] * u[c[i]];
		}
	}

	return ret;
}

double CsrMatrix::mult_vec(size_t irow, const std::vector<double>& u) const{
	const std::vector<size_t>& a = addr();
	const std::vector<size_t>& c = cols();
	const std::vector<double>& v = vals();

	double ret = 0;
	size_t start = a.at(irow);
	size_t end = a.at(irow+1);
	for (size_t i=start; i<end; ++i){
		ret += v[i] * u[c[i]];
	}
	return ret;
}

void CsrMatrix::set_unit_row(size_t irow){
	const std::vector<size_t>& a = addr();
	const std::vector<size_t>& c = cols();

	size_t start = a.at(irow);
	size_t end = a.at(irow+1);
	for (size_t i=start; i<end; ++i){
		_vals[i] = (c[i] == irow) ? 1.0 : 0.0;
	}
}

void CsrMatrix::plus(double k, const CsrMatrix& other){
	if (n_rows() != other.n_rows() || n_nonzeros() != other.n_nonzeros()){
		throw std::runtime_error("can not sum up matrices with different stencils");
	}
	for (size_t i=0; i<n_nonzeros(); ++i){
		_vals[i] += k * other._vals[i];
	}
}


class AmgcMatrixSolver::Impl{
public:
	using param_t = boost::property_tree::ptree;
	using matrix_t = amgcl::backend::crs<double, size_t>;
	using backend_t = amgcl::backend::builtin<double>;
	using solver_t = amgcl::make_solver<
		amgcl::amg<
			backend_t,
			amgcl::runtime::coarsening::wrapper,
			amgcl::runtime::relaxation::wrapper
		>,
		amgcl::runtime::solver::wrapper<backend_t>
	>;

	Impl(matrix_t mat, param_t param): _solver(mat, param), _dim(mat.nrows), _maxit(param.get("solver.maxiter", -1)){ }

	void solve(const std::vector<double>& rhs, std::vector<double>& x) const{
		x.resize(_dim, 0);
		int n_iters;
		double err;
		std::tie(n_iters, err) = _solver(rhs, x);

		if (n_iters >= _maxit){
			if (std::isnan(err)){
				throw std::runtime_error("Sparse matrix solver got NaN");
			}
			std::ostringstream os;
			os << "WARNING: Sparse matrix solution failed to converge with tolerance: "
				<< std::scientific << std::setprecision(2) << err << std::endl;
			std::cout << os.str();
		}
	}
private:
	solver_t _solver;
	int _dim;
	int _maxit;
};


AmgcMatrixSolver::AmgcMatrixSolver(int maxit, double tolerance): _maxit(maxit), _tolerance(tolerance) {}
AmgcMatrixSolver::AmgcMatrixSolver(std::initializer_list<std::pair<std::string, std::string>> amgc_params): AmgcMatrixSolver() {
	for (auto it: amgc_params){
		_params[it.first] = it.second;
	}
}

AmgcMatrixSolver::~AmgcMatrixSolver() = default;

void AmgcMatrixSolver::set_matrix(const CsrMatrix& mat){
	set_matrix(mat, mat.vals());
}

void AmgcMatrixSolver::set_matrix(const CsrStencil& stencil, const std::vector<double>& mat_values){
	Impl::matrix_t amgcl_matrix;
	amgcl_matrix.own_data = false;
	amgcl_matrix.nrows = amgcl_matrix.ncols = stencil.n_rows();
	amgcl_matrix.nnz = stencil.n_nonzeros();
	amgcl_matrix.ptr = const_cast<size_t*>(stencil.addr().data());
	amgcl_matrix.col = const_cast<size_t*>(stencil.cols().data());
	amgcl_matrix.val = const_cast<double*>(mat_values.data());

	Impl::param_t prm;
	prm.put("solver.type", "fgmres");
	prm.put("solver.tol", _tolerance);
	prm.put("solver.maxiter", _maxit);
	prm.put("precond.coarsening.type", "smoothed_aggregation");
	prm.put("precond.relax.type", "spai0");
	//prm.put("precond.relax.type", "gauss_seidel");
	
	for (auto it: _params){
		prm.put(it.first, it.second);
	}

	_pimpl.reset(new Impl(amgcl_matrix, prm));
}

void AmgcMatrixSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
	if (!_pimpl){
		throw std::runtime_error("Matrix was not passed to the solver");
	}
	_pimpl->solve(rhs, ret);
}

void AmgcMatrixSolver::solve_slae(const CsrMatrix& mat, const std::vector<double>& rhs, std::vector<double>& x, int maxit, double eps){
	AmgcMatrixSolver slv(maxit, eps);
	slv.set_matrix(mat);
	slv.solve(rhs, x);
}

std::ostream& operator<<(std::ostream& s, const CsrMatrix& mat){
	int p = 8; 
	for (int i = 0; i < mat.n_rows(); i++)
	{
		for (int j = 0; j < mat.n_rows(); j++)
		{
			if (mat.is_in_stencil(i, j) == false)
			{
				s << std::setw(p) << "*";
			}
			else
			{
				double v = std::round(mat.value(i, j)*10000)/10000.0;
				s << std::setw(p) << v;
			}
		}
		s << std::endl;
	}
	return s;
}

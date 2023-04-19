#include "slae/matrix_solver.hpp"

#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>

using namespace bf;

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
		size_t n_iters;
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
	size_t _dim;
	size_t _maxit;
};


AmgcMatrixSolver::AmgcMatrixSolver(size_t maxit, double tolerance): _maxit(maxit), _tolerance(tolerance) {}
AmgcMatrixSolver::~AmgcMatrixSolver() = default;

void AmgcMatrixSolver::set_matrix(const CsrMatrix& mat){
	set_matrix(mat, mat.vals());
}

void AmgcMatrixSolver::set_matrix(const CsrStencil& stencil, const std::vector<double>& mat_values){
	Impl::matrix_t amgcl_matrix;
	amgcl_matrix.own_data = false;
	amgcl_matrix.nrows = amgcl_matrix.ncols = stencil.n_rows();
	amgcl_matrix.nnz = stencil.n_nonzero();
	amgcl_matrix.ptr = const_cast<size_t*>(stencil.addr().data());
	amgcl_matrix.col = const_cast<size_t*>(stencil.cols().data());
	amgcl_matrix.val = const_cast<double*>(mat_values.data());

	Impl::param_t prm;
	prm.put("solver.type", "fgmres");
	prm.put("solver.tol", _tolerance);
	prm.put("solver.maxiter", _maxit);
	prm.put("precond.coarsening.type", "smoothed_aggregation");
	prm.put("precond.relax.type", "spai0");

	_pimpl.reset(new Impl(amgcl_matrix, prm));
}

void AmgcMatrixSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
	if (!_pimpl){
		throw std::runtime_error("Matrix was not passed to the solver");
	}
	_pimpl->solve(rhs, ret);
}

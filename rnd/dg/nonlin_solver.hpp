#ifndef _DG_NONLIN_SOLVER_HPP_
#define _DG_NONLIN_SOLVER_HPP_

#include "common.hpp"
#include <array>

struct INonlinearSystem2{
	virtual std::array<double, 2> f(double x1, double x2) const = 0;
	virtual std::array<double, 4> jac(double x1, double x2) const = 0;
};

void solve_nonlinear_system(const INonlinearSystem2& sys, double& x1, double& x2, double eps=1e-6, size_t maxit=10'000);


#endif

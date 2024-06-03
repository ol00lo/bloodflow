#include "nonlin_solver.hpp"
#include <sstream>
#include <cmath>

void solve_nonlinear_system(const INonlinearSystem2& sys, double& x1, double& x2, double eps, size_t maxit){

	size_t it = 0;
	double n2 = 0;
	for (it=0; it < maxit; ++it){
		std::array<double, 2> v = sys.f(x1, x2);
		n2 = v[0]*v[0] + v[1]*v[1];
		//std::cout << x1 << " " << x2 << " " << std::sqrt(n2) << std::endl;
		if (n2 < eps*eps){
			break;
		}

		std::array<double, 4> jac= sys.jac(x1, x2);
		double d = jac[0]*jac[3] - jac[1]*jac[2];

		x1 -= (jac[3] * v[0] - jac[1] * v[1])/d;
		x2 -= (-jac[2] * v[0] + jac[0] * v[1])/d;
	}

	if (it >= maxit){
		std::ostringstream oss;
		std::cout << "Warning: ";
		std::cout << "nonlinear system failed to converge in " << maxit << " iterations ";
		std::cout << "till e = " << eps << ". ";
		std::cout << "Norm=" << std::sqrt(n2) << std::endl;
	}
}

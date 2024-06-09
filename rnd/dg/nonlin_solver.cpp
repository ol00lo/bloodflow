#include "nonlin_solver.hpp"
#include <sstream>
#include <cmath>

namespace{
void warning_message(double cur_err, double max_err, size_t maxit){
	std::ostringstream oss;
	std::cout << "Warning: ";
	std::cout << "nonlinear system failed to converge in " << maxit << " iterations ";
	std::cout << "till e = " << max_err << ". ";
	std::cout << "Norm=" << cur_err << std::endl;
}
}

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
		warning_message(std::sqrt(n2), eps, maxit);
	}
}

void solve_nonlinear_system(const INonlinearSystem4& sys, double& x1, double& x2, double& x3, double& x4, double eps, size_t maxit){
	size_t it = 0;
	double n2 = 0;
	for (it=0; it < maxit; ++it){
		std::array<double, 4> v = sys.f(x1, x2, x3, x4);
		n2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
		if (n2 < eps*eps){
			break;
		}

		std::array<double, 16> j = sys.jac(x1, x2, x3, x4);
		double d = j[0]*(j[5]*(j[10]*j[15]-j[11]*j[14])-j[6]*(j[9]*j[15]-j[11]*j[13])+j[7]*(j[9]*j[14]-j[10]*j[13]))-j[1]*(j[4]*(j[10]*j[15]-j[11]*j[14])-j[6]*(j[8]*j[15]-j[11]*j[12])+j[7]*(j[8]*j[14]-j[10]*j[12]))+j[2]*(j[4]*(j[9]*j[15]-j[11]*j[13])-j[5]*(j[8]*j[15]-j[11]*j[12])+j[7]*(j[8]*j[13]-j[9]*j[12]))-j[3]*(j[4]*(j[9]*j[14]-j[10]*j[13])-j[5]*(j[8]*j[14]-j[10]*j[12])+j[6]*(j[8]*j[13]-j[9]*j[12]));
		x1 -= (v[0]*((j[5]*j[10]-j[6]*j[9])*j[15]+(j[7]*j[9]-j[5]*j[11])*j[14]+(j[6]*j[11]-j[7]*j[10])*j[13])+v[1]*((j[2]*j[9]-j[1]*j[10])*j[15]+(j[1]*j[11]-j[3]*j[9])*j[14]+(j[3]*j[10]-j[2]*j[11])*j[13])+v[2]*((j[1]*j[6]-j[2]*j[5])*j[15]+(j[3]*j[5]-j[1]*j[7])*j[14]+(j[2]*j[7]-j[3]*j[6])*j[13])+v[3]*((j[2]*j[5]-j[1]*j[6])*j[11]+(j[1]*j[7]-j[3]*j[5])*j[10]+(j[3]*j[6]-j[2]*j[7])*j[9]))/d;
		x2 -= (v[0]*((j[6]*j[8]-j[4]*j[10])*j[15]+(j[4]*j[11]-j[7]*j[8])*j[14]+(j[7]*j[10]-j[6]*j[11])*j[12])+v[1]*((j[0]*j[10]-j[2]*j[8])*j[15]+(j[3]*j[8]-j[0]*j[11])*j[14]+(j[2]*j[11]-j[3]*j[10])*j[12])+v[2]*((j[2]*j[4]-j[0]*j[6])*j[15]+(j[0]*j[7]-j[3]*j[4])*j[14]+(j[3]*j[6]-j[2]*j[7])*j[12])+v[3]*((j[0]*j[6]-j[2]*j[4])*j[11]+(j[3]*j[4]-j[0]*j[7])*j[10]+(j[2]*j[7]-j[3]*j[6])*j[8]))/d;
		x3 -= (v[0]*((j[4]*j[9]-j[5]*j[8])*j[15]+(j[7]*j[8]-j[4]*j[11])*j[13]+(j[5]*j[11]-j[7]*j[9])*j[12])+v[1]*((j[1]*j[8]-j[0]*j[9])*j[15]+(j[0]*j[11]-j[3]*j[8])*j[13]+(j[3]*j[9]-j[1]*j[11])*j[12])+v[2]*((j[0]*j[5]-j[1]*j[4])*j[15]+(j[3]*j[4]-j[0]*j[7])*j[13]+(j[1]*j[7]-j[3]*j[5])*j[12])+v[3]*((j[1]*j[4]-j[0]*j[5])*j[11]+(j[0]*j[7]-j[3]*j[4])*j[9]+(j[3]*j[5]-j[1]*j[7])*j[8]))/d;
		x4 -= (v[0]*((j[5]*j[8]-j[4]*j[9])*j[14]+(j[4]*j[10]-j[6]*j[8])*j[13]+(j[6]*j[9]-j[5]*j[10])*j[12])+v[1]*((j[0]*j[9]-j[1]*j[8])*j[14]+(j[2]*j[8]-j[0]*j[10])*j[13]+(j[1]*j[10]-j[2]*j[9])*j[12])+v[2]*((j[1]*j[4]-j[0]*j[5])*j[14]+(j[0]*j[6]-j[2]*j[4])*j[13]+(j[2]*j[5]-j[1]*j[6])*j[12])+v[3]*((j[0]*j[5]-j[1]*j[4])*j[10]+(j[2]*j[4]-j[0]*j[6])*j[9]+(j[1]*j[6]-j[2]*j[5])*j[8]))/d;
	}

	if (it >= maxit){
		warning_message(std::sqrt(n2), eps, maxit);
	}
}

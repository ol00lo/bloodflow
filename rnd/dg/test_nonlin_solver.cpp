#include "common_test.hpp"
#include "nonlin_solver.hpp"

namespace{

struct Test1_NonlinSystem: public INonlinearSystem2{
	std::array<double, 2> f(double x1, double x2) const override{
		return {
			x1*x1 + x2*x2 - 25,
			x1*x1 - x2*x2 + 7
		};
	}

	std::array<double, 4> jac(double x1, double x2) const override{
		return {
			2*x1, 2*x2,
			2*x1, -2*x2
		};
	}
};

}

TEST_CASE("Two quadratic equations system", "[nonlin-solver][1]"){
	Test1_NonlinSystem sys;

	{
		double x1 = 2;
		double x2 = 2;
		solve_nonlinear_system(sys, x1, x2);
		CHECK(x1 == Approx(3).margin(1e-6));
		CHECK(x2 == Approx(4).margin(1e-6));
	}

	{
		double x1 = -2;
		double x2 = -2;
		solve_nonlinear_system(sys, x1, x2);
		CHECK(x1 == Approx(-3).margin(1e-6));
		CHECK(x2 == Approx(-4).margin(1e-6));
	}
}

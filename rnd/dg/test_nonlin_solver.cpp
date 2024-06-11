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


namespace{

struct Test2_NonlinSystem: public INonlinearSystem4{
	std::array<double, 4> f(double x1, double x2, double x3, double x4) const override{
		return {
			x1*x1 + x2*x2 + x3 + x3*x4 - 4,
			x1*x1 - x2*x2 - x3*x4 + 1,
			x1*x4 - 1,
			x2*x3 - 1,
		};
	}

	std::array<double, 16> jac(double x1, double x2, double x3, double x4) const override{
		return {
			2*x1, 2*x2, 1+x4, x3,
			2*x1, -2*x2, -x4, -x3,
			x4, 0, 0, x1,
			0, x3, x2, 0
		};
	}
};

}

TEST_CASE("Four quadratic equations system", "[nonlin-solver][2]"){
	Test2_NonlinSystem sys;

	{
		double x1 = 1.1;
		double x2 = 1.1;
		double x3 = 0.9;
		double x4 = 0.9;
		solve_nonlinear_system(sys, x1, x2, x3, x4);
		CHECK(x1 == Approx(1).margin(1e-6));
		CHECK(x2 == Approx(1).margin(1e-6));
		CHECK(x3 == Approx(1).margin(1e-6));
		CHECK(x4 == Approx(1).margin(1e-6));
	}
}


namespace{

struct Test3_NonlinSystem: public INonlinearSystem6{
	std::array<double, 6> f(double x1, double x2, double x3, double x4, double x5, double x6) const override{
		return {
			x1*x1 + x2*x2 + x3 + x3*x4 - x1*x5 + 2*x6 - 5,
			x1*x1 - x2*x2 - x3*x4 + x3*x5 - x6 + 1,
			x1*x4 + x5*x6 - 2,
			x2*x3 - x5*x6,
			x1 + x6 - 2,
			x2 - x5,
		};
	}

	std::array<double, 36> jac(double x1, double x2, double x3, double x4, double x5, double x6) const override{
		return {
			2*x1-x5, 2*x2, 1+x4, x3, -x1, 2,
			2*x1, -2*x2, -x4+x5, -x3, x3, -1,
			x4, 0, 0, x1, x6, x5,
			0, x3, x2, 0, -x6, -x5,
			1, 0, 0, 0, 0, 1,
			0, 1, 0, 0, -1, 0
		};
	}
};

}

TEST_CASE("Six quadratic equations system", "[nonlin-solver][3]"){
	Test3_NonlinSystem sys;

	{
		double x1 = 1.1;
		double x2 = 0.9;
		double x3 = 1.1;
		double x4 = 0.9;
		double x5 = 1.1;
		double x6 = 0.9;
		solve_nonlinear_system(sys, x1, x2, x3, x4, x5, x6);
		CHECK(x1 == Approx(1).margin(1e-6));
		CHECK(x2 == Approx(1).margin(1e-6));
		CHECK(x3 == Approx(1).margin(1e-6));
		CHECK(x4 == Approx(1).margin(1e-6));
		CHECK(x5 == Approx(1).margin(1e-6));
		CHECK(x6 == Approx(1).margin(1e-6));
	}
}

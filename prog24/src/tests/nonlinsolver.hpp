#ifndef _DG_NONLIN_SOLVER_HPP_
#define _DG_NONLIN_SOLVER_HPP_
#include <array>

struct INonlinearSystem2
{
    virtual std::array<double, 2> f(double x1, double x2) const = 0;
    virtual std::array<double, 4> jac(double x1, double x2) const = 0;
};
struct INonlinearSystem4
{
    virtual std::array<double, 4> f(double x1, double x2, double x3, double x4) const = 0;
    virtual std::array<double, 16> jac(double x1, double x2, double x3, double x4) const = 0;
};
struct INonlinearSystem6
{
    virtual std::array<double, 6> f(double x1, double x2, double x3, double x4, double x5, double x6) const = 0;
    virtual std::array<double, 36> jac(double x1, double x2, double x3, double x4, double x5, double x6) const = 0;
};

void solve_nonlinear_system(const INonlinearSystem2& sys, double& x1, double& x2, double eps = 1e-12,
                            size_t maxit = 10'000);
void solve_nonlinear_system(const INonlinearSystem4& sys, double& x1, double& x2, double& x3, double& x4,
                            double eps = 1e-12, size_t maxit = 10'000);
void solve_nonlinear_system(const INonlinearSystem6& sys, double& x1, double& x2, double& x3, double& x4, double& x5,
                            double& x6, double eps = 1e-12, size_t maxit = 10'000);

#endif

#include "tests/nonlinsolver.hpp"
#include <iostream>
#include <cmath>
#include <sstream>

void solve_nonlinear_system(const INonlinearSystem2& sys, double& x1, double& x2, double eps, size_t maxit)
{

    size_t it = 0;
    double n2 = 0;
    for (it = 0; it < maxit; ++it)
    {
        std::array<double, 2> v = sys.f(x1, x2);
        n2 = v[0] * v[0] + v[1] * v[1];
        // std::cout << x1 << " " << x2 << " " << std::sqrt(n2) << std::endl;
        if (n2 < eps * eps)
        {
            break;
        }

        std::array<double, 4> jac = sys.jac(x1, x2);
        double d = jac[0] * jac[3] - jac[1] * jac[2];

        x1 -= (jac[3] * v[0] - jac[1] * v[1]) / d;
        x2 -= (-jac[2] * v[0] + jac[0] * v[1]) / d;
    }

    if (it >= maxit)
    {
        std::ostringstream oss;
        std::cout << "Warning: ";
        std::cout << "nonlinear system failed to converge in " << maxit << " iterations ";
        std::cout << "till e = " << eps << ". ";
        std::cout << "Norm=" << std::sqrt(n2) << std::endl;
    }
}
double find_d_3_3(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8)
{
    return a0 * (a4 * a8 - a7 * a5) + a6 * (a1 * a5 - a4 * a2) + a3 * (a7 * a2 - a1 * a8);
}
std::array<double, 16> inverse4_4(std::array<double, 16> a, double det)
{
    std::array<double, 16> res;
    double a11 = a[0];
    double a12 = a[1];
    double a13 = a[2];
    double a14 = a[3];
    double a21 = a[4];
    double a22 = a[5];
    double a23 = a[6];
    double a24 = a[7];
    double a31 = a[8];
    double a32 = a[9];
    double a33 = a[10];
    double a34 = a[11];
    double a41 = a[12];
    double a42 = a[13];
    double a43 = a[14];
    double a44 = a[15];
    res[0] = +find_d_3_3(a22, a23, a24, a32, a33, a34, a42, a43, a44) / det;
    res[1] = -find_d_3_3(a21, a23, a24, a31, a33, a34, a41, a43, a44) / det;
    res[2] = +find_d_3_3(a21, a22, a24, a31, a32, a34, a41, a42, a44) / det;
    res[3] = -find_d_3_3(a21, a22, a23, a31, a32, a33, a41, a42, a43) / det;

    res[4] = -find_d_3_3(a12, a13, a14, a32, a33, a34, a42, a43, a44) / det;
    res[5] = +find_d_3_3(a11, a13, a14, a31, a33, a34, a41, a43, a44) / det;
    res[6] = -find_d_3_3(a11, a12, a14, a31, a32, a34, a41, a42, a44) / det;
    res[7] = +find_d_3_3(a11, a12, a13, a31, a32, a33, a41, a42, a43) / det;

    res[8] = +find_d_3_3(a12, a13, a14, a22, a23, a24, a42, a43, a44) / det;
    res[9] = -find_d_3_3(a11, a13, a14, a21, a23, a24, a41, a43, a44) / det;
    res[10] = +find_d_3_3(a11, a12, a14, a21, a22, a24, a41, a42, a44) / det;
    res[11] = -find_d_3_3(a11, a12, a13, a21, a22, a23, a41, a42, a43) / det;

    res[12] = -find_d_3_3(a12, a13, a14, a22, a23, a24, a32, a33, a34) / det;
    res[13] = +find_d_3_3(a11, a13, a14, a21, a23, a24, a31, a33, a34) / det;
    res[14] = -find_d_3_3(a11, a12, a14, a21, a22, a24, a31, a32, a34) / det;
    res[15] = +find_d_3_3(a11, a12, a13, a21, a22, a23, a31, a32, a33) / det;
    return res;
}

void solve_nonlinear_system(const INonlinearSystem4& sys, double& x1, double& x2, double& x3, double& x4, double eps,
                            size_t maxit)
{
    size_t it = 0;
    double n2 = 0;
    for (it = 0; it < maxit; ++it)
    {
        std::array<double, 4> v = sys.f(x1, x2, x3, x4);
        n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
        std::cout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << std::sqrt(n2) << std::endl;
        if (n2 < eps * eps)
        {
            break;
        }
        std::array<double, 16> jac = sys.jac(x1, x2, x3, x4);
        double d = jac[0] * find_d_3_3(jac[5], jac[6], jac[7], jac[9], jac[10], jac[11], jac[13], jac[14], jac[15]) +
                   jac[4] * find_d_3_3(jac[1], jac[2], jac[3], jac[9], jac[10], jac[11], jac[13], jac[14], jac[15]) +
                   jac[8] * find_d_3_3(jac[1], jac[2], jac[3], jac[5], jac[6], jac[7], jac[13], jac[14], jac[15]) +
                   jac[12] * find_d_3_3(jac[1], jac[2], jac[3], jac[5], jac[6], jac[7], jac[9], jac[10], jac[11]);
        std::array<double, 16> jac_inv = inverse4_4(jac, d);

        x1 -= (jac_inv[0] * v[0] + jac_inv[1] * v[1] + jac_inv[2] * v[2] + jac_inv[3] * v[3]);
        x2 -= (jac_inv[4] * v[0] + jac_inv[5] * v[1] + jac_inv[6] * v[2] + jac_inv[7] * v[3]);
        x3 -= (jac_inv[8] * v[0] + jac_inv[9] * v[1] + jac_inv[10] * v[2] + jac_inv[11] * v[3]);
        x4 -= (jac_inv[12] * v[0] + jac_inv[13] * v[1] + jac_inv[14] * v[2] + jac_inv[15] * v[3]);
    }

    if (it >= maxit)
    {
        std::ostringstream oss;
        std::cout << "Warning: ";
        std::cout << "nonlinear system failed to converge in " << maxit << " iterations ";
        std::cout << "till e = " << eps << ". ";
        std::cout << "Norm=" << std::sqrt(n2) << std::endl;
    }
}
void solve_nonlinear_system(const INonlinearSystem6& sys, double& x1, double& x2, double& x3, double& x4, double& x5,
                            double& x6, double eps, size_t maxit)
{
    size_t it = 0;
    double n2 = 0;
    for (it = 0; it < maxit; ++it)
    {
    }

    if (it >= maxit)
    {
        std::ostringstream oss;
        std::cout << "Warning: ";
        std::cout << "nonlinear system failed to converge in " << maxit << " iterations ";
        std::cout << "till e = " << eps << ". ";
        std::cout << "Norm=" << std::sqrt(n2) << std::endl;
    }
}

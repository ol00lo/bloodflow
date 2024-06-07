#include "tests/nonlinsolver.hpp"
#include <iostream>
#include <cmath>
#include <vector>
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
    return a0 * (a4 * a8 - a7 * a5) + a1 * (a6 * a5 - a3 * a8) + a2 * (a7 * a3 - a4 * a6);
}

double determinant(std::vector<std::vector<double>> matrix, int size)
{
    int det = 0;
    int sign = 1;
    for (int i = 0; i < size; i++)
    {
        std::vector<std::vector<double>> cofactor(size - 1);
        for (int j = 0; j < size - 1; j++)
            cofactor[j].resize(size-1);

        int sub_i = 0, sub_j = 0;
        for (int j = 1; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                if (k == i)
                {
                    continue;
                }
                cofactor[sub_i][sub_j] = matrix[j][k];
                sub_j++;
            }
            sub_i++;
            sub_j = 0;
        }
        det += sign * matrix[0][i] * determinant(cofactor, size - 1);
        sign = -sign;
    }
    return det;
}

std::array<double, 16> inverse4_4(std::array<double, 16> a, double det)
{
    std::array<double, 16> res;
    res[0] = +find_d_3_3(a[5], a[6], a[7], a[9], a[10], a[11], a[13], a[14], a[15]) / det;
    res[1] = -find_d_3_3(a[1], a[2], a[3], a[9], a[10], a[11], a[13], a[14], a[15]) / det;
    res[2] = +find_d_3_3(a[1], a[2], a[3], a[5], a[6], a[7], a[13], a[14], a[15]) / det;
    res[3] = -find_d_3_3(a[1], a[2], a[3], a[5], a[6], a[7], a[9], a[10], a[11]) / det;

    res[4] = -find_d_3_3(a[4], a[6], a[7], a[8], a[10], a[11], a[12], a[14], a[15]) / det;
    res[5] = +find_d_3_3(a[0], a[2], a[3], a[8], a[10], a[11], a[12], a[14], a[15]) / det;
    res[6] = -find_d_3_3(a[0], a[2], a[3], a[4], a[6], a[7], a[12], a[14], a[15]) / det;
    res[7] = +find_d_3_3(a[0], a[2], a[3], a[4], a[6], a[7], a[8], a[10], a[11]) / det;

    res[8] = +find_d_3_3(a[4], a[5], a[7], a[8], a[9], a[11], a[12], a[13], a[15]) / det;
    res[9] = -find_d_3_3(a[0], a[1], a[3], a[8], a[9], a[11], a[12], a[13], a[15]) / det;
    res[10] = +find_d_3_3(a[0], a[1], a[3], a[4], a[5], a[7], a[12], a[13], a[15]) / det;
    res[11] = -find_d_3_3(a[0], a[1], a[3], a[4], a[5], a[7], a[8], a[9], a[11]) / det;

    res[12] = -find_d_3_3(a[4], a[5], a[6], a[8], a[9], a[10], a[12], a[13], a[14]) / det;
    res[13] = +find_d_3_3(a[0], a[1], a[2], a[8], a[9], a[10], a[12], a[13], a[14]) / det;
    res[14] = -find_d_3_3(a[0], a[1], a[2], a[4], a[5], a[6], a[12], a[13], a[14]) / det;
    res[15] = +find_d_3_3(a[0], a[1], a[2], a[4], a[5], a[6], a[8], a[9], a[10]) / det;
    return res;
}

std::array<double, 36> inverse6_6(std::vector<std::vector<double>> a, double det)
{
    std::array<double, 36> res;
    res[0] = determinant({a[11], a[12], a[13], a[14], a[15]}, 4);
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
        std::array<double, 6> v = sys.f(x1, x2, x3, x4, x5, x6);
        n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]+ v[4]*v[4]+v[5]*v[5];
        std::cout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << std::sqrt(n2) << std::endl;
        if (n2 < eps * eps)
        {
            break;
        }
        std::array<double, 36> jac = sys.jac(x1, x2, x3, x4, x5, x6);
        std::vector<std::vector<double>> matrix = {{jac[0], jac[1], jac[2], jac[3], jac[4], jac[5]},
                                                   {jac[6], jac[7], jac[8], jac[9], jac[10], jac[11]},
                                                   {jac[12], jac[13], jac[14], jac[15], jac[16], jac[17]},
                                                   {jac[18], jac[19], jac[20], jac[21], jac[22], jac[23]},
                                                   {jac[24], jac[25], jac[26], jac[27], jac[28], jac[29]},
                                                   {jac[30], jac[31], jac[32], jac[33], jac[34], jac[35]}};
        double d = determinant(matrix, 0.01);
        std::array<double, 36> jac_inv = inverse6_6(matrix, d);

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

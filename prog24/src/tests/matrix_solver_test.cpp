#include "bflow/matrix_solver.hpp"
#include "catch.hpp"
#include <iostream>
#include <sstream>

using namespace bflow;

TEST_CASE("test solve matrix 1", "[solve_matrix1]")
{
    std::vector<int> addr;
    std::vector<int> cols;
    std::vector<double> vals;
    std::vector<double> rhs;
    addr = {0, 1, 2, 3};
    cols = {0, 1, 2};
    vals = {1, 1, 1};
    rhs = {2, 4, 5};
    std::vector<double> u;
    CsrMatrix mat(addr, cols, vals);
    AmgcMatrixSolver a;
    a.set_matrix(mat);
    a.solve(rhs, u);
    CHECK(u[0] == 2);
    CHECK(u[1] == 4);
    CHECK(u[2] == 5);
    addr = {0, 1, 2, 3};
    cols = {1, 0, 0};
    vals = {1, 1, 1};
    rhs = {1, 1, 1};
    CsrMatrix mat1(addr, cols, vals);
    CHECK_THROWS(a.set_matrix(mat1));
    cols = {0, 1, 2};
    addr = {0, 1, 1, 3};
    CsrMatrix mat2(addr, cols, vals);
    CHECK_THROWS(a.set_matrix(mat2));

    addr = {0, 2, 3};
    cols = {0, 1, 1};
    vals = {2, 1, 3};
    rhs = {7.0 / 3, 1};
    CsrMatrix mat3(addr, cols, vals);
    a.set_matrix(mat3);
    a.solve(rhs, u);
    CHECK(u[0] == 1);
    CHECK(u[1] == Approx(1.0 / 3));
}
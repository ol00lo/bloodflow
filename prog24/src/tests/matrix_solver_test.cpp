#include "bflow/matrix_solver.hpp"
#include "catch.hpp"
#include <sstream>
#include <iostream>

using namespace bflow;

TEST_CASE("test solve matrix 1", "[solve_matrix1]")
{
    std::vector<int> addr = {0,2,3};
    std::vector<int> cols = {0,1,1};
    std::vector<double> vals = {2,1,3};
    std::vector<double> rhs = {7.0/3,1};
    std::vector<double> u;
    CsrMatrix mat(addr, cols, vals);

    AmgcMatrixSolver a;
    a.set_matrix(mat);
    a.solve(rhs, u);
    for (auto s : u)
    {
        std::cout << s << std::endl;
    }
}
#include "bflow/fem_grid.hpp"
#include "bflow/graph_grid.hpp"
#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/vessel_graph.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include <fstream>

#include <iomanip>
using namespace bflow;

namespace
{
double exact(double x, double t = 0)
{
    constexpr double eps = 1e-6;
    if (t > 0)
        return exact(x - t);
    else
        return (x > eps && x < 0.8 - eps) ? 1.0 : 0.0;
}

double exact1(double x, double t = 0)
{
    constexpr double pi = 3.1415926535;
    return sin(2 * pi * (x - t));
}

double norm2(FemGrid grid, std::vector<double> u, double time)
{
    auto lv = grid.load_vector();
    double I = 0;
    double gamma = 0;
    for (size_t i = 0; i < grid.n_points(); ++i)
    {
        double r = u[i] - exact(grid.node(i), time);
        I += lv[i] * (r * r);
        gamma += grid.h(0);
    }
    return std::sqrt(I / gamma);
}

void print_matrix_full(const CsrMatrix& mat, std::ostream& s = std::cout)
{
    for (size_t i = 0; i < mat.n_rows(); i++)
    {
        for (size_t j = 0; j < mat.n_rows(); j++)
        {
            if (mat.is_in_stencil(i, j) == false)
            {
                s << std::setw(6) << "*";
            }
            else
            {
                s << std::setw(6) << round(mat.value(i, j) * 1000) / 1000;
            }
        }
        s << std::endl;
    }
}
} // namespace

TEST_CASE("Transport equation, upwind", "[upwind-transport]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {3.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1, 2);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1);
    double tau = grid.h(0) / 2;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.block_transport_matrix();
    // upwind coupling
    for (size_t ielem = 0; ielem < grid.n_elements(); ielem++)
    {
        std::vector<int> lg = grid.tab_elem_nodes(ielem);
        if (ielem > 0)
        {
            // left
            size_t iaddr = transport.find_index(lg[0], lg[0] - 1);
            transport.vals()[iaddr] += 1;
        }
        {
            // right
            size_t iaddr = transport.find_index(lg[1], lg[1]);
            transport.vals()[iaddr] -= 1;
        }
    }
    CsrMatrix lhs = mass;
    CsrMatrix rhs_mat = mass;
    double theta = 0.5; // crank-nicolson
    // double theta = 1.5;  // adams-bashforth
    for (size_t i = 0; i < mass.n_nonzeros(); ++i)
    {
        lhs.vals()[i] -= theta * tau * transport.vals()[i];
        rhs_mat.vals()[i] += (1 - theta) * tau * transport.vals()[i];
    }
    // left boundary condition
    lhs.set_unit_row(0);

    // matrix solver
    AmgcMatrixSolver slv;
    slv.set_matrix(lhs);

    // initial conditions
    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact1(grid.node(i));
    }

    NonstatGridSaver saver(grid1, nodes_coo, "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(u, "data");

    for (double t = tau; t <= 2.0; t += tau)
    {
        // assemble rhs
        std::vector<double> rhs = rhs_mat.mult_vec(u);
        // left boundary condition
        rhs[0] = exact1(0, t);

        slv.solve(rhs, u);

        // std::cout << t << "  ";
        // std::cout << norm2(grid, u, t) << std::endl;

        saver.new_time_step(t);
        saver.save_vtk_point_data(u, "data");
    }
    // CHECK(u[50] == Approx(1.071884924).margin(1e-6));
    CHECK(norm2(grid, u, 2.0) == Approx(0.284369755).margin(1e-4));
}

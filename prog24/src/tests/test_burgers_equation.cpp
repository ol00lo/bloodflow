#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/fem_grid.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include <fstream>
using namespace bflow;
namespace
{
double exact(double x, double t = 0)
{
    return x / (t + 1);
}
} // namespace

TEST_CASE("Inviscid Burgers equation, explicit", "[Burgers-inviscid-explicit]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1);
    double tau = grid.h() / 2;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.block_transport_matrix();
    // upwind coupling
    for (int ielem = 0; ielem < grid.n_elements(); ielem++)
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
    // left boundary condition
    lhs.set_unit_row(0);

    AmgcMatrixSolver slv;
    slv.set_matrix(lhs);

    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact(grid.node(i));
    }
    NonstatGridSaver saver(grid1, nodes_coo, "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(u, "data");

    double time = 0;
    while (time + tau < 2.0 + 1e-6)
    {
        time += tau;
        std::cout << time << std::endl;
        // assemble rhs
        std::vector<double> flux(grid.n_nodes());
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            flux[i] = u[i] * u[i] / 2;
        }
        std::vector<double> fa = transport.mult_vec(flux);
        std::vector<double> rhs = mass.mult_vec(u);
        for (size_t i = 0; i < rhs.size(); ++i)
        {
            rhs[i] += tau * fa[i];
        }

        // left boundary condition
        rhs[0] = 0.0;

        slv.solve(rhs, u);

        saver.new_time_step(time);
        saver.save_vtk_point_data(u, "data");
    }
    CHECK(u.back() == Approx(0.3274573644).margin(1e-6));
}

TEST_CASE("Inviscid Burgers equation, implicit", "[Burgers-inviscid-implicit][amg]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1);
    double tau = grid.h() / 2;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.block_transport_matrix();
    // upwind coupling
    for (int ielem = 0; ielem < grid.n_elements(); ielem++)
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
    std::vector<double> load_vector = mass.mult_vec(std::vector<double>(grid.n_nodes(), 1.0));

    AmgcMatrixSolver slv(10'000, 1e-14);

    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact(grid.node(i));
    }
    NonstatGridSaver saver(grid1, nodes_coo, "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(u, "data");
    size_t maxit = 1000;
    double maxeps = 1e-14;
    double time = 0;
    while (time + tau < 2.0 + 1e-6)
    {
        time += tau;
        std::cout << time << std::endl;

        // assemble rhs
        std::vector<double> rhs = mass.mult_vec(u);
        rhs[0] = 0.0;

        for (size_t it = 0; it < maxit; ++it)
        {
            // assemble lhs
            CsrMatrix lhs = mass;
            for (size_t i = 0; i < lhs.n_rows(); ++i)
            {
                for (size_t a = lhs.addr()[i]; a < lhs.addr()[i + 1]; ++a)
                {
                    size_t col = lhs.cols()[a];
                    lhs.vals()[a] -= tau * u[col] / 2 * transport.vals()[a];
                }
            }
            // left boundary condition
            lhs.set_unit_row(0);

            // compute residual
            if (it > 0)
            {
                std::vector<double> r = lhs.mult_vec(u);
                double err = 0;
                for (size_t i = 0; i < r.size(); ++i)
                {
                    double diff = lhs.mult_vec(i, u) - rhs[i];
                    err += diff * diff * load_vector[i];
                }
                err = std::sqrt(err);
                if (err < maxeps)
                {
                    std::cout << "converged in " << it << " iterations" << std::endl;
                    break;
                }
                if (it == maxit - 1)
                {
                    std::cout << "Failed to converge with err=" << err << std::endl;
                }
            }

            slv.set_matrix(lhs);
            slv.solve(rhs, u);
        }

        saver.new_time_step(time);
        saver.save_vtk_point_data(u, "data");
    }
    CHECK(u.back() == Approx(0.3396867099).margin(1e-6));
}

TEST_CASE("Inviscid Burgers equation, implicit, iterflux", "[Burgers-inviscid-implicit-iterflux][amg]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1);
    double tau = grid.h() / 5;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix block_transport = grid.block_transport_matrix();
    CsrMatrix coupled_transport = grid.coupling_transport_matrix();
    std::vector<double> load_vector = mass.mult_vec(std::vector<double>(grid.n_nodes(), 1.0));

    AmgcMatrixSolver slv(10'000, 1e-14);

    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact(grid.node(i));
    }
    NonstatGridSaver saver(grid1, nodes_coo, "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(u, "data");

    size_t maxit = 10000;
    double maxeps = 1e-14;
    double time = 0;
    while (time + tau < 2.0 + 1e-6)
    {
        time += tau;
        std::cout << time << std::endl;

        // assemble rhs
        std::vector<double> rhs1 = mass.mult_vec(u);

        for (size_t it = 0; it < maxit; ++it)
        {
            // assemble lhs
            CsrMatrix lhs = mass;
            for (size_t i = 0; i < lhs.n_rows(); ++i)
            {
                for (size_t a = lhs.addr()[i]; a < lhs.addr()[i + 1]; ++a)
                {
                    size_t col = lhs.cols()[a];
                    lhs.vals()[a] -= tau * u[col] / 2 * block_transport.vals()[a];
                }
            }

            // assemble rhs
            std::vector<double> rhs(rhs1);
            std::vector<double> u2(grid.n_nodes());
            for (size_t i = 0; i < u2.size(); ++i)
                u2[i] = u[i] * u[i];

            //// var1: goes to rhs
            // for (size_t i=0; i<grid.n_nodes(); ++i){
            //        rhs[i] -= tau * coupled_transport.mult_vec(i, u2) / 2.0;
            //}

            // var2: goes to the lhs diagonal
            for (size_t i = 0; i < lhs.n_rows(); ++i)
            {
                double denum = (u[i] == 0) ? 1e-16 : u[i];
                double v = tau * coupled_transport.mult_vec(i, u2) / 2.0 / denum;
                for (size_t a = lhs.addr()[i]; a < lhs.addr()[i + 1]; ++a)
                {
                    size_t col = lhs.cols()[a];
                    if (col == i)
                    {
                        lhs.vals()[a] += v;
                        break;
                    }
                }
            }

            //// var3: goes to the both sides
            // for (size_t i=0; i<lhs.n_rows(); ++i){
            //        double v = coupled_transport.mult_vec(i, u2) / 2.0;
            //        size_t il, ir;
            //        if (i == 0 || i == lhs.n_rows() - 1){
            //                il = ir = i;
            //        } else if (i%2 == 0){
            //                il = i - 1;
            //                ir = i;
            //        } else {
            //                il = i;
            //                ir = i+1;
            //        }
            //        for (size_t a=lhs.addr()[i]; a < lhs.addr()[i+1]; ++a){
            //                size_t col = lhs.cols()[a];
            //                if (col == ir){
            //                        double d = (u[col] == 0) ? 1e-16 : u[col];
            //                        lhs.vals()[a] += tau * v / 2 / d;
            //                }
            //                if (col == il){
            //                        double d = (u[col] == 0) ? 1e-16 : u[col];
            //                        lhs.vals()[a] += tau * v / 2 / d;
            //                }
            //        }
            //}

            // left boundary condition
            lhs.set_unit_row(0);
            rhs[0] = 0;

            // compute residual
            if (it > 0)
            {
                std::vector<double> r = lhs.mult_vec(u);
                double err = 0;
                for (size_t i = 0; i < r.size(); ++i)
                {
                    double diff = lhs.mult_vec(i, u) - rhs[i];
                    err += diff * diff * load_vector[i];
                }
                err = std::sqrt(err);
                if (err < maxeps)
                {
                    std::cout << "converged in " << it << " iterations" << std::endl;
                    break;
                }
                if (it == maxit - 1)
                {
                    std::cout << "Failed to converge with err=" << err << std::endl;
                }
            }

            slv.set_matrix(lhs);
            slv.solve(rhs, u);
        }

        saver.new_time_step(time);
        saver.save_vtk_point_data(u, "data");
    }
    CHECK(u.back() == Approx(0.3360909317).margin(1e-6));
}

TEST_CASE("Inviscid Burgers equation, cn", "[Burgers-inviscid-cn][amg]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1);
    double tau = grid.h() / 2;
    double theta = 0.5;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.block_transport_matrix();
    // upwind coupling
    for (int ielem = 0; ielem < grid.n_elements(); ielem++)
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
    std::vector<double> load_vector = mass.mult_vec(std::vector<double>(grid.n_nodes(), 1.0));

    AmgcMatrixSolver slv(10'000, 1e-14);

    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact(grid.node(i));
    }
    NonstatGridSaver saver(grid1, nodes_coo, "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(u, "data");

    size_t maxit = 1000;
    double maxeps = 1e-14;
    double time = 0;
    while (time + tau < 2.0 + 1e-6)
    {
        time += tau;
        std::cout << time << std::endl;

        // assemble rhs
        std::vector<double> rhs = mass.mult_vec(u);
        std::vector<double> flux(grid.n_nodes());
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            flux[i] = u[i] * u[i] / 2;
        }
        flux = transport.mult_vec(flux);
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            rhs[i] += (1 - theta) * tau * flux[i];
        }
        rhs[0] = 0.0;

        for (size_t it = 0; it < maxit; ++it)
        {
            // assemble lhs
            CsrMatrix lhs = mass;
            for (size_t i = 0; i < lhs.n_rows(); ++i)
            {
                for (size_t a = lhs.addr()[i]; a < lhs.addr()[i + 1]; ++a)
                {
                    size_t col = lhs.cols()[a];
                    lhs.vals()[a] -= theta * tau * u[col] / 2 * transport.vals()[a];
                }
            }
            // left boundary condition
            lhs.set_unit_row(0);

            // compute residual
            if (it > 0)
            {
                std::vector<double> r = lhs.mult_vec(u);
                double err = 0;
                for (size_t i = 0; i < r.size(); ++i)
                {
                    double diff = lhs.mult_vec(i, u) - rhs[i];
                    err += diff * diff * load_vector[i];
                }
                err = std::sqrt(err);
                if (err < maxeps)
                {
                    std::cout << "converged in " << it << " iterations" << std::endl;
                    break;
                }
                if (it == maxit - 1)
                {
                    std::cout << "Failed to converge with err=" << err << std::endl;
                }
            }

            slv.set_matrix(lhs);
            slv.solve(rhs, u);
        }

        saver.new_time_step(time);
        saver.save_vtk_point_data(u, "data");
    }
    CHECK(u.back() == Approx(0.3335733261).margin(1e-6));
}

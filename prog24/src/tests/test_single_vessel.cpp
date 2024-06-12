#include "bflow/assembler_flux.hpp"
#include "bflow/fem_grid.hpp"
#include "bflow/flux_calculator.hpp"
#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include <fstream>
#include <iomanip>
using namespace bflow;

namespace
{
double q_inflow1(double t)
{
    return 1e-6 * exp(-1e4 * (t - 0.05) * (t - 0.05));
};
}
TEST_CASE("Single vessel, inviscid2, ver1", "[single-vessel-inviscid-explicit1]")
{
    ProblemData data;
    double time = 0;
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {data.L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, data.L/100, 1);
    FemGrid grid(grid1);
    double tau = grid.h() / 100;
    std::vector<ElementBoundaryFluxes> upwind_fluxes(grid.n_elements());

    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    NonstatGridSaver saver(grid1, nodes_coo, "bebe");

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes(), data.area0);
    std::vector<double> pressure(grid.n_nodes(), 0.0);
    saver.new_time_step(0);
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(velocity, "velocity");
    saver.save_vtk_point_data(pressure, "pressure");

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
    upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, [&time]() { return time; }, 0));
    for (size_t i = 1; i < grid.n_points() - 1; ++i)
    {
        upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i - 1, i));
    }
    upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements() - 1));

    // prepare matrix solver
    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix tran = grid.block_transport_matrix();
    AmgcMatrixSolver slv;
    slv.set_matrix(mass);
    double t = 0;
    while (time < 0.2 - 1e-6)
    {
        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  Q=" << data.q_inflow(time) << std::endl;

        // nodewise fluxes
        std::vector<double> flux_a(grid.n_nodes());
        std::vector<double> flux_u(grid.n_nodes());
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            flux_a[i] = data.flux_a(area[i], velocity[i]);
            flux_u[i] = data.flux_u(area[i], velocity[i]);
        }

        // upwind fluxes
        for (auto c : upwind_flux_calculator)
        {
            c->compute(area, velocity, upwind_fluxes);
        }

        // assemble rhs
        // E
        std::vector<double> rhs_a = mass.mult_vec(area);
        std::vector<double> rhs_u = mass.mult_vec(velocity);
        // +tau*T
        std::vector<double> tran_a = tran.mult_vec(flux_a);
        std::vector<double> tran_u = tran.mult_vec(flux_u);
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            rhs_a[i] -= tau * tran_a[i];
            rhs_u[i] -= tau * tran_u[i];
        }
        // + coupling
        for (size_t ielem = 0; ielem < grid.n_elements(); ++ielem)
        {
            size_t node0 = grid.tab_elem_nodes(ielem)[0];
            size_t node1 = grid.tab_elem_nodes(ielem)[1];
            rhs_a[node0] += tau * upwind_fluxes[ielem].a_x0;
            rhs_a[node1] -= tau * upwind_fluxes[ielem].a_x1;
            rhs_u[node0] += tau * upwind_fluxes[ielem].u_x0;
            rhs_u[node1] -= tau * upwind_fluxes[ielem].u_x1;
        }
        // solve
        slv.solve(rhs_a, area);
        slv.solve(rhs_u, velocity);

        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            pressure[i] = data.pressure(area[i]);
        }
        t += tau;
        if (t > 0.005-1e-6)
        {
            t = 0;
            saver.new_time_step(time);
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(velocity, "velocity");
            saver.save_vtk_point_data(pressure, "pressure");
        }
    }

    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(24.3958).margin(1e-3));
}

TEST_CASE("Single vessel, inviscid, implicit", "[single-vessel-inviscid-implicit]")
{
    ProblemData data;
    data.L = 1;
    data.recompute();
    double time = 0;
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {data.L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, data.L / 30, 1);
    FemGrid grid(grid1);
    double tau = grid.h() / 100;

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes(), data.area0);
    std::vector<double> pressure(grid.n_nodes(), 0.0);
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        pressure[i] = data.pressure(area[i]);
    }
    std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
    upwind_flux_calculator[0].reset(new InflowQFluxCalculator(
        grid, data, [&time]() { return q_inflow1(time); }, 0));
    for (size_t i = 1; i < grid.n_points() - 1; ++i)
    {
        upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i - 1, i));
    }
    upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements() - 1));

    AssemblerFlux assem(grid, &data, upwind_flux_calculator);
    AmgcMatrixSolver slv;

    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    NonstatGridSaver saver(grid1, nodes_coo, "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(velocity, "velocity");
    saver.save_vtk_point_data(pressure, "pressure");

    size_t iter_max = 100;
    double t = 0.0;
    while (time < 0.2 - 1e-6)
    {
        // assemble right hand side
        std::vector<double> rhs1_a = assem.mass().mult_vec(area);
        std::vector<double> rhs1_u = assem.mass().mult_vec(velocity);

        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  Q=" << data.q_inflow(time) << std::endl;

        double err_a = 1e6;
        double err_u = 1e6;

        for (size_t it = 0; it < iter_max; ++it)
        {
            // 1. --- Area equation
            assem.actualize_fluxes(time, area, velocity);
            // 1.1 LHS
            CsrMatrix lhs_a = assem.mass();
            CsrMatrix block_u_transport = assem.block_u_transport();
            for (size_t i = 0; i < lhs_a.n_nonzeros(); ++i)
            {
                lhs_a.vals()[i] -= tau * block_u_transport.vals()[i];
            }
            // 1.2 RHS
            std::vector<double> rhs_a(rhs1_a);
            // coupling flux
            std::vector<double> coupling_flux_ua = assem.coupling_flux_ua();
            for (size_t i = 0; i < rhs_a.size(); ++i)
            {
                rhs_a[i] -= tau * coupling_flux_ua[i];
            }
            // 1.3 residual
            if (it > 0)
            {
                err_a = assem.compute_residual(lhs_a, rhs_a, area);
            }
            std::cout << std::endl;
            // 1.4 Solve
            slv.set_matrix(lhs_a);
            slv.solve(rhs_a, area);
            // 2. ---- Velocity equation
            assem.actualize_fluxes(time, area, velocity);

            // 2.1 LHS
            CsrMatrix lhs_u = assem.mass();
            block_u_transport = assem.block_u_transport();
            for (size_t i = 0; i < lhs_u.n_nonzeros(); ++i)
            {
                lhs_u.vals()[i] -= 0.5 * tau * block_u_transport.vals()[i];
            }

            // 2.2 RHS
            std::vector<double> rhs_u(rhs1_u);
            // coupling flux
            std::vector<double> coupling_flux_u2 = assem.coupling_flux_u2();
            std::vector<double> coupling_flux_p = assem.coupling_flux_p();
            for (size_t i = 0; i < rhs_a.size(); ++i)
            {
                rhs_u[i] -= tau * (0.5 * coupling_flux_u2[i] + coupling_flux_p[i] / data.rho);
            }
            // block p
            std::vector<double> p(grid.n_nodes());
            for (size_t i = 0; i < grid.n_nodes(); ++i)
            {
                p[i] = data.pressure(area[i]);
            }
            auto tmp = assem.block_transport().mult_vec(p);
            for (size_t i = 0; i < rhs_a.size(); ++i)
            {
                rhs_u[i] -= tau * tmp[i] / data.rho;
            }

            // 2.3 residual

            if (it > 0)
            {
                err_u = assem.compute_residual(lhs_u, rhs_u, velocity);
            }

            // 2.4 Solve
            slv.set_matrix(lhs_u);

            slv.solve(rhs_u, velocity);
            // break conditions
            if (it > 0)
            {
                std::cout << err_u << " " << err_a << std::endl;
                if (std::max(err_u, err_a) < 1e-12)
                {
                    std::cout << "converged in " << it << " iterations" << std::endl;
                    break;
                }
                else if (it == iter_max - 1)
                {
                    std::cout << "Warning: internal iterations did not converge: " << err_u << ", " << err_a
                              << std::endl;
                }
            }
        }

        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            pressure[i] = data.pressure(area[i]);
        }
        t += tau;
        if (t>0.005-1e-6)
        {
            t = 0;
            saver.new_time_step(time);
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(velocity, "velocity");
            saver.save_vtk_point_data(pressure, "pressure");
        }
    }

    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(14.0287775452).margin(1e-3));
}

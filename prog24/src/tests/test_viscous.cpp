#include "bflow/flux_calculator.hpp"
#include "bflow/assembler_flux.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/fem_grid.hpp"
#include "bflow/macros.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"

namespace
{
double q_inflow1(double t)
{
    return 1e-6 * exp(-1e4 * (t - 0.05) * (t - 0.05));
};
} // namespace

TEST_CASE("Single vessel, viscous", "[single-vessel-viscous]")
{
    ProblemData data;
    data.mu = 0.04;
    data.recompute();

    double L = 10;
    double time = 0;
    double theta = 0.5;

    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, L / 1000, 1);
    FemGrid grid(grid1);
    double tau = grid.h() / 50;

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
    NonstatGridSaver saver(grid1, generate_nodes_coo(gr1), "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(velocity, "velocity");
    saver.save_vtk_point_data(pressure, "pressure");

    size_t iter_max = 10;

    CsrMatrix block_u_transport;
    std::vector<double> coupling_flux_ua;
    std::vector<double> coupling_flux_u2;
    std::vector<double> coupling_flux_p;
    std::vector<double> tmp;

    double t = 0;
    while (time < 1.5 - 1e-6)
    {
        // assemble right hand side
        assem.actualize_fluxes(time, area, velocity);
        block_u_transport = assem.block_u_transport();
        // 1. ---- Area equation
        std::vector<double> rhs1_a = assem.mass().mult_vec(area);
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            rhs1_a[i] += (1 - theta) * tau * block_u_transport.mult_vec(i, area);
        }
        coupling_flux_ua = assem.coupling_flux_ua();
        for (size_t i = 0; i < rhs1_a.size(); ++i)
        {
            rhs1_a[i] -= (1 - theta) * tau * coupling_flux_ua[i];
        }
        // 2. ---- Velocity equation
        std::vector<double> rhs1_u = assem.mass().mult_vec(velocity);
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] += 0.5 * (1 - theta) * tau * block_u_transport.mult_vec(i, velocity);
        }
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            pressure[i] = data.pressure(area[i]);
        }
        tmp = assem.block_transport().mult_vec(pressure);
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] += (1 - theta) * tau * tmp[i] / data.rho;
        }
        coupling_flux_u2 = assem.coupling_flux_u2();
        coupling_flux_p = assem.coupling_flux_p();
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] -= (1 - theta) * tau * (coupling_flux_u2[i] / 2 + coupling_flux_p[i] / data.rho);
        }
        tmp = assem.viscous_matrix().mult_vec(velocity);
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            rhs1_u[i] += (1 - theta) * tau * tmp[i];
        }

        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  Q=" << q_inflow1(time) << std::endl;

        double err_a = 1e6;
        double err_u = 1e6;
        double err_max = 1e-10;

        for (size_t it = 0; it < iter_max; ++it)
        {
            // 1. --- Area equation
            assem.actualize_fluxes(time, area, velocity);
            // 1.1 LHS
            CsrMatrix lhs_a = assem.mass();
            lhs_a.plus(-theta * tau, assem.block_u_transport());
            // 1.2 RHS
            std::vector<double> rhs_a(rhs1_a);
            // coupling flux
            coupling_flux_ua = assem.coupling_flux_ua();
            for (size_t i = 0; i < rhs_a.size(); ++i)
            {
                rhs_a[i] -= theta * tau * coupling_flux_ua[i];
            }
            // 1.3 residual
            if (it > 0)
            {
                err_a = assem.compute_residual(lhs_a, rhs_a, area);
            }

            // 1.4 Solve
            slv.set_matrix(lhs_a);
            slv.solve(rhs_a, area);

            // 2. ---- Velocity equation
            assem.actualize_fluxes(time, area, velocity);
            // 2.1 LHS
            CsrMatrix lhs_u = assem.mass();
            lhs_u.plus(-0.5 * theta * tau, assem.block_u_transport());
            lhs_u.plus(-theta * tau, assem.viscous_matrix());
            // 2.2 RHS
            std::vector<double> rhs_u(rhs1_u);
            // coupling flux
            coupling_flux_u2 = assem.coupling_flux_u2();
            coupling_flux_p = assem.coupling_flux_p();
            for (size_t i = 0; i < rhs_u.size(); ++i)
            {
                rhs_u[i] -= theta * tau * (0.5 * coupling_flux_u2[i] + coupling_flux_p[i] / data.rho);
            }
            // block p
            for (size_t i = 0; i < grid.n_nodes(); ++i)
            {
                pressure[i] = data.pressure(area[i]);
            }
            tmp = assem.block_transport().mult_vec(pressure);
            for (size_t i = 0; i < rhs_u.size(); ++i)
            {
                rhs_u[i] += theta * tau * tmp[i] / data.rho;
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
                // std::cout << err_u << " " << err_a << std::endl;
                if (std::max(err_u, err_a) < err_max)
                {
                    // std::cout << "converged in " << it << " iterations" << std::endl;
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
        if (t > 0.05 - 1e-6)
        {
            t = 0;
            saver.new_time_step(time);
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(velocity, "velocity");
            saver.save_vtk_point_data(pressure, "pressure");
        }

    }
    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(15.2062197472).margin(1e-3));
}

#include "bflow/assembler_flux.hpp"
#include "bflow/debug/printer.hpp"
#include "bflow/fem_grid.hpp"
#include "bflow/flux_calculator.hpp"
#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include <fstream>

TEST_CASE("Bifurcated vessel", "[bifurcated-vessel]")
{
    ProblemData data1;
    data1.area0 = ProblemData::pi * 1e-4 / 4;
    data1.rho = 1000;
    data1.h = 1;
    data1.E = 3.0 * data1.area0 * 324970 / 4.0 / std::sqrt(ProblemData::pi);
    data1.recompute();

    ProblemData data2;
    data2.area0 = ProblemData::pi * 1e-4 / 6.0 / 4.0;
    data2.rho = 1000;
    data2.h = 1;
    data2.E = 3.0 * data1.area0 * 796020 / 4.0 / std::sqrt(ProblemData::pi);
    data2.recompute();
    data1.report();
    data2.report();
    throw;

    double time = 0;
    double theta = 0.5;
    double L = 0.2;
    // size_t k3 = 50;
    size_t k3 = 10;

    std::vector<std::vector<int>> node = {{0}, {0, 1, 2}, {1}, {2}};
    std::vector<double> ed = {L, L, L};
    // std::vector<std::vector<int>> node = {{0}, {0, 1}, {1}};
    // std::vector<double> ed = {L, L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, L / k3, 1);
    FemGrid grid(grid1);
    std::vector<int> cell_types(grid.n_elements(), 2);
    for (size_t i = 0; i < k3; ++i)
    {
        cell_types[i] = 1;
    }
    double tau = L / k3 / 20;

    size_t monitoring_node_A = grid.closest_node(0, 0.0 * L);
    size_t monitoring_node_B = grid.closest_node(0, 0.5 * L);
    size_t monitoring_node_C = grid.closest_node(0, 0.9999 * L);
    size_t monitoring_node_D = grid.closest_node(1, 0.5 * L);
    size_t monitoring_node_E = grid.closest_node(1, 0.9999 * L);
    std::vector<double> monitor_A, monitor_B, monitor_C, monitor_D, monitor_E;

    std::vector<const ProblemData*> data(grid.n_elements());
    for (size_t icell = 0; icell < grid.n_elements(); ++icell)
    {
        data[icell] = (cell_types[icell] == 1) ? &data1 : &data2;
    }

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> flux_calculators(grid.n_points());
    // inflow
    auto getau = [&time]() -> std::pair<double, double> {
        constexpr double time0 = 0.05;
        constexpr double C = 5000;
        constexpr double a = 0.0001 * ProblemData::pi / 4;

        double t = time - time0;
        double u = 0.01 * exp(-C * t * t);

        return {a, u};
    };
    flux_calculators[0].reset(new InflowAUFluxCalculator(grid, data1, getau, 0));
    // internal
    for (size_t i = 1; i < k3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data1, i - 1, i));
    for (size_t i = k3 + 1; i < 2 * k3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i - 1, i));
    for (size_t i = 2 * k3 + 1; i < 3 * k3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i - 1, i));
    // outflow
    flux_calculators[2 * k3].reset(new OutflowFluxCalculator(grid, data2, 2 * k3 - 1));
    flux_calculators[3 * k3].reset(new OutflowFluxCalculator(grid, data2, 3 * k3 - 1));
    // bifurcation
    flux_calculators[k3].reset(new Bifurcation3FluxCalculator(grid, data1, data2, data2, k3 - 1, k3, 2 * k3));

    // assembler
    AssemblerFlux assem(grid, data, flux_calculators);

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        area[i] = (cell_types[grid.tab_node_elem(i)] == 1) ? data1.area0 : data2.area0;
    }

    AmgcMatrixSolver slv;

    NonstatGridSaver saver(grid1, generate_nodes_coo(gr1), "bifurcated-vessel");
    saver.new_time_step(0);
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(velocity, "velocity");

    size_t iter_max = 10;

    CsrMatrix block_u_transport;
    std::vector<double> coupling_flux_ua;
    std::vector<double> coupling_flux_u2;
    std::vector<double> coupling_flux_p;
    std::vector<double> tmp;
    double t = 0;
    // while (time < 0.5 - 1e-12)
    while (time < 0.05 - 1e-12)
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
        // dbg::print_matrix_full(assem.block_transport());
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
        tmp = assem.block_transport().mult_vec(assem.pressure());
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] += (1 - theta) * tau * tmp[i] / data1.rho;
        }
        coupling_flux_u2 = assem.coupling_flux_u2();
        coupling_flux_p = assem.coupling_flux_p();
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] -= (1 - theta) * tau * (coupling_flux_u2[i] / 2 + coupling_flux_p[i] / data1.rho);
        }
        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  U=" << getau().second << std::endl;

        double err_a = 1e6;
        double err_u = 1e6;
        double err_max = 1e-8;

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
            // 2.2 RHS
            std::vector<double> rhs_u(rhs1_u);
            // coupling flux
            coupling_flux_u2 = assem.coupling_flux_u2();
            coupling_flux_p = assem.coupling_flux_p();
            for (size_t i = 0; i < rhs_u.size(); ++i)
            {
                rhs_u[i] -= theta * tau * (0.5 * coupling_flux_u2[i] + coupling_flux_p[i] / data1.rho);
            }
            // block p
            tmp = assem.block_transport().mult_vec(assem.pressure());
            for (size_t i = 0; i < rhs_u.size(); ++i)
            {
                rhs_u[i] += theta * tau * tmp[i] / data1.rho;
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
        t += tau;
        if (t > 0.005 - 1e-6)
        {
            t = 0;
            saver.new_time_step(time);
            // saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(velocity, "velocity");
            saver.save_vtk_point_data(assem.pressure(area), "pressure");
            saver.save_vtk_point_data(assem.w1(), "w1");
            saver.save_vtk_point_data(assem.w2(), "w2");
            saver.save_vtk_point_data(assem.normalized_area(area), "area");
        }

        std::vector<double> pressure = assem.pressure(area);

        monitor_A.push_back(pressure[monitoring_node_A]);
        monitor_B.push_back(pressure[monitoring_node_B]);
        monitor_C.push_back(pressure[monitoring_node_C]);
        monitor_D.push_back(pressure[monitoring_node_D]);
        monitor_E.push_back(pressure[monitoring_node_E]);
    }

    std::vector<double> pressure = assem.pressure(area);
    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(6.262212285).margin(1e-3));

    time_value_vtk(tau, monitor_A, "monitor_A.vtk");
    time_value_vtk(tau, monitor_B, "monitor_B.vtk");
    time_value_vtk(tau, monitor_C, "monitor_C.vtk");
    time_value_vtk(tau, monitor_D, "monitor_D.vtk");
    time_value_vtk(tau, monitor_E, "monitor_E.vtk");
}

TEST_CASE("Bifurcated vessel1", "[bifurcated-vessel1]")
{
    ProblemData data1;
    data1.area0 = ProblemData::pi * 5e-3 * 5e-3;
    data1.rho = 1050;
    data1.h = 0.15e-3;
    data1.E = 1e6;
    data1.recompute();

    ProblemData data2 = data1;
    data2.area0 = ProblemData::pi * 4e-3 * 4e-3;
    data2.recompute();

    ProblemData data2_alter = data2;
    // data2_alter.area0 = data2.area0 / 2;
    // data2_alter.E = 10*data2.E;
    data2_alter.recompute();

    ProblemData data3 = data1;
    data3.area0 = ProblemData::pi * 3e-3 * 3e-3;
    data3.recompute();

    data1.report();
    data2.report();
    data2_alter.report();
    data3.report();

    double time = 0;
    double theta = 0.5;
    double L1 = 1;
    double L2 = 0.8;
    double L3 = 0.5;
    // size_t n1 = 50;
    size_t n1 = 10;
    size_t n2 = L2 / (L1 / n1);
    size_t n3 = L3 / (L1 / n1);
    double tau = L1 / n1 / 200;
    double save_tau = 0.01;
    size_t iter_max = 10;

    std::vector<std::vector<int>> node = {{0}, {0, 1, 2}, {1, 3, 4}, {2, 5, 6}, {3}, {4}, {5}, {6}};
    std::vector<double> ed = {L1, L2, L2, L3, L3, L3, L3};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, L1 / n1, 1);
    if (grid1.n_elem() != n1 + 2 * n2 + 4 * n3)
    {
        throw std::runtime_error("invalid grid partition");
    }
    FemGrid grid(grid1);
    std::vector<int> cell_types(grid.n_elements());
    for (size_t i = 0; i < grid.n_elements(); ++i)
    {
        if (i < n1)
            cell_types[i] = 1;
        else if (i < n1 + n2)
            cell_types[i] = 2;
        else if (i < n1 + 2 * n2)
            cell_types[i] = 20;
        else
            cell_types[i] = 3;
    }

    std::vector<const ProblemData*> data(grid.n_elements());
    for (size_t icell = 0; icell < grid.n_elements(); ++icell)
    {
        switch (cell_types[icell])
        {
        case 1:
            data[icell] = &data1;
            break;
        case 2:
            data[icell] = &data2;
            break;
        case 20:
            data[icell] = &data2_alter;
            break;
        case 3:
            data[icell] = &data3;
            break;
        default:
            throw std::runtime_error("ERROR");
        }
    }

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> flux_calculators(grid.n_points());
    // inflow
    auto getau = [&time, &data1]() -> std::pair<double, double> {
        double bpm = 60.0;
        double T = 60.0 / bpm;
        double t0 = T / 2;
        double q0 = 2e-5;
        double tm = time;
        while (tm > T)
            tm -= T;
        double Q = (tm < t0) ? q0 * sin(ProblemData::pi * tm / t0) : 0;
        double a = data1.area0;
        double u = Q / a;
        return {a, u};
    };
    flux_calculators[0].reset(new InflowAUFluxCalculator(grid, data1, getau, 0));
    // internal
    for (size_t i = 1; i < n1; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data1, i - 1, i));

    for (size_t i = n1 + 1; i < n1 + n2; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i - 1, i));

    for (size_t i = n1 + n2 + 1; i < n1 + 2 * n2; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data2_alter, i - 1, i));

    for (size_t i = n1 + 2 * n2 + 1; i < n1 + 2 * n2 + n3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data3, i - 1, i));

    for (size_t i = n1 + 2 * n2 + n3 + 1; i < n1 + 2 * n2 + 2 * n3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data3, i - 1, i));

    for (size_t i = n1 + 2 * n2 + 2 * n3 + 1; i < n1 + 2 * n2 + 3 * n3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data3, i - 1, i));

    for (size_t i = n1 + 2 * n2 + 3 * n3 + 1; i < n1 + 2 * n2 + 4 * n3; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data3, i - 1, i));

    // outflow
    flux_calculators[n1 + 2 * n2 + n3].reset(new OutflowFluxCalculator(grid, data3, n1 + 2 * n2 + n3 - 1));
    flux_calculators[n1 + 2 * n2 + 2 * n3].reset(new OutflowFluxCalculator(grid, data3, n1 + 2 * n2 + 2 * n3 - 1));
    flux_calculators[n1 + 2 * n2 + 3 * n3].reset(new OutflowFluxCalculator(grid, data3, n1 + 2 * n2 + 3 * n3 - 1));
    flux_calculators[n1 + 2 * n2 + 4 * n3].reset(new OutflowFluxCalculator(grid, data3, n1 + 2 * n2 + 4 * n3 - 1));

    // bifurcation
    flux_calculators[n1].reset(new Bifurcation3FluxCalculator(grid, data1, data2, data2_alter, n1 - 1, n1, n1 + n2));
    flux_calculators[n1 + n2].reset(
        new Bifurcation3FluxCalculator(grid, data2, data3, data3, n1 + n2 - 1, n1 + 2 * n2, n1 + 2 * n2 + n3));
    flux_calculators[n1 + 2 * n2].reset(new Bifurcation3FluxCalculator(grid, data2_alter, data3, data3, n1 + 2 * n2 - 1,
                                                                       n1 + 2 * n2 + 2 * n3, n1 + 2 * n2 + 3 * n3));

    // assembler
    AssemblerFlux assem(grid, data, flux_calculators);

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        switch (cell_types[grid.tab_node_elem(i)])
        {
        case 1:
            area[i] = data1.area0;
            break;
        case 2:
            area[i] = data2.area0;
            break;
        case 20:
            area[i] = data2_alter.area0;
            break;
        case 3:
            area[i] = data3.area0;
            break;
        default:
            throw std::runtime_error("ERROR");
        }
    }

    AmgcMatrixSolver slv;

    NonstatGridSaver saver(grid1, generate_nodes_coo(gr1), "bifurcated-vessel");
    auto save_data = [&]() {
        assem.actualize_fluxes(time, area, velocity);
        saver.new_time_step(time);
        saver.save_vtk_point_data(velocity, "velocity");
        saver.save_vtk_point_data(assem.pressure(), "pressure");
        saver.save_vtk_point_data(assem.w1(), "w1");
        saver.save_vtk_point_data(assem.w2(), "w2");
        saver.save_vtk_point_data(assem.normalized_area(), "area");
    };

    save_data();

    size_t q_monitor_node_11 = grid.tab_point_nodes(n1 + 2 * n2 + 1 * n3)[0];
    size_t q_monitor_node_12 = grid.tab_point_nodes(n1 + 2 * n2 + 2 * n3)[0];
    size_t q_monitor_node_21 = grid.tab_point_nodes(n1 + 2 * n2 + 3 * n3)[0];
    size_t q_monitor_node_22 = grid.tab_point_nodes(n1 + 2 * n2 + 4 * n3)[0];
    std::vector<double> q_monitor_1;
    std::vector<double> q_monitor_2;

    CsrMatrix block_u_transport;
    std::vector<double> coupling_flux_ua;
    std::vector<double> coupling_flux_u2;
    std::vector<double> coupling_flux_p;
    std::vector<double> tmp;

    double t = 0;
    // while (time < 2.0 - 1e-12)
    while (time < 0.5 - 1e-12)
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
        // dbg::print_matrix_full(assem.block_transport());
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
        tmp = assem.block_transport().mult_vec(assem.pressure());
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] += (1 - theta) * tau * tmp[i] / data1.rho;
        }
        coupling_flux_u2 = assem.coupling_flux_u2();
        coupling_flux_p = assem.coupling_flux_p();
        for (size_t i = 0; i < rhs1_u.size(); ++i)
        {
            rhs1_u[i] -= (1 - theta) * tau * (coupling_flux_u2[i] / 2 + coupling_flux_p[i] / data1.rho);
        }
        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  U=" << getau().second << std::endl;

        double err_a = 1e6;
        double err_u = 1e6;
        double err_max = 1e-8;

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
            // 2.2 RHS
            std::vector<double> rhs_u(rhs1_u);
            // coupling flux
            coupling_flux_u2 = assem.coupling_flux_u2();
            coupling_flux_p = assem.coupling_flux_p();
            for (size_t i = 0; i < rhs_u.size(); ++i)
            {
                rhs_u[i] -= theta * tau * (0.5 * coupling_flux_u2[i] + coupling_flux_p[i] / data1.rho);
            }
            // block p
            tmp = assem.block_transport().mult_vec(assem.pressure());
            for (size_t i = 0; i < rhs_u.size(); ++i)
            {
                rhs_u[i] += theta * tau * tmp[i] / data1.rho;
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
        t += tau;
        if (t > save_tau - 1e-6)
        {
            t = 0;
            save_data();
        }

        double q11 = area[q_monitor_node_11] * velocity[q_monitor_node_11];
        double q12 = area[q_monitor_node_12] * velocity[q_monitor_node_12];
        double q21 = area[q_monitor_node_21] * velocity[q_monitor_node_21];
        double q22 = area[q_monitor_node_22] * velocity[q_monitor_node_22];
        q_monitor_1.push_back(q11 + q12);
        q_monitor_2.push_back(q21 + q22);
    }

    std::vector<double> pressure = assem.pressure(area);
    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(549.4302681774).margin(1e-3));
    CHECK(q_monitor_1.back() == Approx(0.0000006744).margin(1e-9));

    time_value_vtk(tau, q_monitor_1, "monitor_q1.vtk");
    time_value_vtk(tau, q_monitor_2, "monitor_q2.vtk");
}

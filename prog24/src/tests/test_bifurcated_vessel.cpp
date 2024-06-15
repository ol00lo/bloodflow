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
namespace
{
void time_value_vtk(double tau, const std::vector<double>& v, std::string filename)
{
    std::ofstream fs(filename);
    size_t n_nodes = v.size();
    size_t n_cells = v.size() - 1;
    fs << "# vtk DataFile Version 3.0" << std::endl;
    fs << "DG" << std::endl;
    fs << "ASCII" << std::endl;
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
    fs << "POINTS " << n_nodes << " double" << std::endl;
    for (size_t i = 0; i < n_nodes; ++i)
    {
        fs << i * tau << " 0 0" << std::endl;
    }

    // Cells
    fs << "CELLS  " << n_cells << "   " << 3 * n_cells << std::endl;
    for (size_t ielem = 0; ielem < n_cells; ++ielem)
    {
        fs << 2 << " " << ielem << " " << ielem + 1 << std::endl;
    }
    fs << "CELL_TYPES  " << n_cells << std::endl;
    for (size_t i = 0; i < n_cells; ++i)
        fs << 3 << std::endl;

    // Data
    fs << "POINT_DATA " << n_nodes << std::endl;
    fs << "SCALARS area  double 1" << std::endl;
    fs << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < n_nodes; ++i)
    {
        fs << v[i] << std::endl;
    }
    fs.close();
}
void PRVK(std::vector<double> v)
{
    std::cout << std::endl;
    for (auto x : v)
        std::cout << x << " ";
}
} // namespace

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
            saver.save_vtk_point_data(assem.area(area), "area");
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
    ProblemData data3;
    data3.area0 = ProblemData::pi * 2e-3 * 2e-3;
    data3.rho = 1050;
    data3.h = 0.15e-3;
    data3.E = 6e6;
    data3.recompute();

    ProblemData data2;
    data2.area0 = ProblemData::pi * 3e-3 * 3e-3;
    data2.rho = 1050;
    data2.h = 0.15e-3;
    data2.E = 6e6;
    data2.recompute();

    ProblemData data1;
    data1.area0 = ProblemData::pi * 4e-3 * 4e-3;
    data1.rho = 1050;
    data1.h = 0.15e-3;
    data1.E = 6e6;
    data1.recompute();

    double time = 0;
    double theta = 0.5;
    double L1 = 0.2;
    double L2 = 0.15;
    double L3 = 0.1;
    size_t n1 = 20;
    size_t n2 = 15;
    size_t n3 = 10;

    std::vector<std::vector<int>> node = {{0}, {0, 1, 2}, {1, 3, 4}, {2, 5, 6}, {3}, {4}, {5}, {6}};
    std::vector<double> ed = {L1, L2, L2, L3, L3, L3, L3};
    // std::vector<std::vector<int>> node = {{0}, {0, 1}, {1}};
    // std::vector<double> ed = {L, L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, L1/n1, 1);
    FemGrid grid(grid1);
    std::vector<int> cell_types(grid.n_elements(), 3);
    for (size_t i = 0; i < n1+n2; ++i)
    {
        cell_types[i] = 2;
    }
    for (size_t i = 0; i < n1; ++i)
    {
        cell_types[i] = 1;
    }
    double tau = L1 / n1 / 100;

    // size_t monitoring_node_A = grid.closest_node(0, 0.0 * L1);
    // size_t monitoring_node_B = grid.closest_node(0, 0.5 * L2);
    // size_t monitoring_node_C = grid.closest_node(0, 0.9999 * L2);
    // size_t monitoring_node_D = grid.closest_node(1, 0.5 * L2);
    // size_t monitoring_node_E = grid.closest_node(1, 0.9999 * L2);
    // std::vector<double> monitor_A, monitor_B, monitor_C, monitor_D, monitor_E;

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
        double q0 = 2e-5;
        double Q = (time < 0.5) ? q0 * sin(2 * ProblemData::pi * time) : 0;
        double a = data1.area0;
        double  u = Q / a;
        return {a, u};
    };
    flux_calculators[0].reset(new InflowAUFluxCalculator(grid, data1, getau, 0));
    // internal
    for (size_t i = 1; i < n1; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data1, i - 1, i));

    for (size_t i = n1 + 1; i < n1 + n2; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i - 1, i));

    for (size_t i = n1 + n2 + 1; i < n1 + 2 * n2; ++i)
        flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i - 1, i));

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
    flux_calculators[n1].reset(new Bifurcation3FluxCalculator(grid, data1, data2, data2, n1 - 1, n1, n1 + n2));
    flux_calculators[n1 + n2].reset(
        new Bifurcation3FluxCalculator(grid, data2, data3, data3, n1 + n2 - 1, n1 + 2 * n2, n1 + 2 * n2 + n3));
    flux_calculators[n1 + 2 * n2].reset(new Bifurcation3FluxCalculator(grid, data2, data3, data3, n1 + 2 * n2 - 1,
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
        case 3:
            area[i] = data3.area0;
            break;
        default:
            throw std::runtime_error("ERROR");
        }
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
    // 
    while (time < 2.0 - 1e-12)
    //while (time < 0.05 - 1e-12)
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
            saver.save_vtk_point_data(assem.area(area), "area");
        }

        std::vector<double> pressure = assem.pressure(area);

        // monitor_A.push_back(pressure[monitoring_node_A]);
        // monitor_B.push_back(pressure[monitoring_node_B]);
        // monitor_C.push_back(pressure[monitoring_node_C]);
        // monitor_D.push_back(pressure[monitoring_node_D]);
        // monitor_E.push_back(pressure[monitoring_node_E]);
    }

    std::vector<double> pressure = assem.pressure(area);
    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(6.262212285).margin(1e-3));

    // time_value_vtk(tau, monitor_A, "monitor_A.vtk");
    // time_value_vtk(tau, monitor_B, "monitor_B.vtk");
    // time_value_vtk(tau, monitor_C, "monitor_C.vtk");
    // time_value_vtk(tau, monitor_D, "monitor_D.vtk");
    // time_value_vtk(tau, monitor_E, "monitor_E.vtk");
}

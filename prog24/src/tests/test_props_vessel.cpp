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
double p_inflow2(double t)
{
    constexpr double T = 0.33;
    // constexpr double T = 0.01;
    return (T / 2 - t > 0) ? 2000 * sin(2 * ProblemData::pi * t / T) : 0;
};
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
} // namespace
TEST_CASE("Single vessel, different beta properties", "[2props-vessel]")
{
    ProblemData data1;
    data1.area0 = ProblemData::pi / 4.0;
    data1.rho = 1;
    data1.h = 1;
    data1.E = 84'628.5 * std::sqrt(ProblemData::pi);
    data1.recompute();

    ProblemData data2;
    data2.area0 = ProblemData::pi / 4.0;
    data2.rho = 1;
    data2.h = 1;
    data2.E = data1.E * 100;
    data2.recompute();

    double time = 0;
    double theta = 0.5;
    double L = 15;
    size_t k3 = 10;

    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, L / 3 / k3, 1);
    FemGrid grid(grid1);
    std::vector<int> cell_types(grid.n_elements(), 1);
    for (size_t i = k3; i < 2 * k3; ++i)
    {
        cell_types[i] = 2;
    }
    double tau = grid.h() / 50000;

    size_t monitoring_node1 = grid.closest_node(0.25 * L);
    size_t monitoring_node2 = grid.closest_node(0.5 * L);
    size_t monitoring_node3 = grid.closest_node(0.75 * L);
    std::vector<double> monitor1, monitor2, monitor3;

    std::vector<const ProblemData*> data(grid.n_elements());
    for (size_t icell = 0; icell < grid.n_elements(); ++icell)
    {
        data[icell] = (cell_types[icell] == 1) ? &data1 : &data2;
    }

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> flux_calculators(grid.n_points());
    flux_calculators[0].reset(new InflowPFluxCalculator_NonReflecting(
        grid, data1, [&time]() { return p_inflow2(time); }, 0));
    for (size_t i = 1; i < grid.n_points() - 1; ++i)
    {
        size_t cell_left = i - 1;
        size_t cell_right = i;
        ProblemData* data_left = (cell_types[cell_left] == 1) ? &data1 : &data2;
        ProblemData* data_right = (cell_types[cell_right] == 1) ? &data1 : &data2;
        if (data_left == data_right)
        {
            flux_calculators[i].reset(new InternalFluxCalculator(grid, *data_left, cell_left, cell_right));
        }
        else
        {
            flux_calculators[i].reset(
                new MergingFluxCalculator(grid, *data_left, *data_right, cell_left, cell_right, 1e-8));
        }
    }
    flux_calculators.back().reset(new OutflowFluxCalculator(grid, data1, grid.n_elements() - 1));
    AssemblerFlux assem(grid, data, flux_calculators);

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes(), data1.area0);

    AmgcMatrixSolver slv;
    NonstatGridSaver saver(grid1, generate_nodes_coo(gr1), "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(velocity, "velocity");

    size_t iter_max = 1;

    CsrMatrix block_u_transport;
    std::vector<double> coupling_flux_ua;
    std::vector<double> coupling_flux_u2;
    std::vector<double> coupling_flux_p;
    std::vector<double> tmp;
    double t = 0;
    // while (time < 0.25 - 1e-12){
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
        std::cout << "  P=" << p_inflow2(time) << std::endl;

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
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(velocity, "velocity");
        }

        std::vector<double> pressure = assem.pressure(area);
        monitor1.push_back(pressure[monitoring_node1]);
        monitor2.push_back(pressure[monitoring_node2]);
        monitor3.push_back(pressure[monitoring_node3]);
    }

    std::vector<double> pressure = assem.pressure(area);
    CHECK(pressure[20] == Approx(1521.4710994706).margin(1e-3));

    time_value_vtk(tau, monitor1, "monitor1.vtk");
    time_value_vtk(tau, monitor2, "monitor2.vtk");
    time_value_vtk(tau, monitor3, "monitor3.vtk");
}

TEST_CASE("Single vessel, different beta properties, power==2", "[2props-vessel2]")
{
    ProblemData data1;
    data1.area0 = ProblemData::pi / 4.0;
    data1.rho = 1;
    data1.h = 1;
    data1.E = 84'628.5 * std::sqrt(ProblemData::pi);
    data1.recompute();

    ProblemData data2;
    data2.area0 = ProblemData::pi / 4.0;
    data2.rho = 1;
    data2.h = 1;
    data2.E = data1.E * 100;
    data2.recompute();

    double time = 0;
    double theta = 0.5;
    double L = 15;
    size_t k3 = 10;

    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, L / 3 / k3, 2);
    FemGrid grid(grid1);
    std::vector<int> cell_types(grid.n_elements(), 1);
    for (size_t i = k3; i < 2 * k3; ++i)
    {
        cell_types[i] = 2;
    }
    double tau = grid.h() / 1000000;

    size_t monitoring_node1 = grid.closest_node(0.25 * L);
    size_t monitoring_node2 = grid.closest_node(0.5 * L);
    size_t monitoring_node3 = grid.closest_node(0.75 * L);
    std::vector<double> monitor1, monitor2, monitor3;

    std::vector<const ProblemData*> data(grid.n_elements());
    for (size_t icell = 0; icell < grid.n_elements(); ++icell)
    {
        data[icell] = (cell_types[icell] == 1) ? &data1 : &data2;
    }

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> flux_calculators(grid.n_points());
    flux_calculators[0].reset(new InflowPFluxCalculator_NonReflecting(
        grid, data1, [&time]() { return p_inflow2(time); }, 0));
    for (size_t i = 1; i < grid.n_points() - 1; ++i)
    {
        size_t cell_left = i - 1;
        size_t cell_right = i;
        ProblemData* data_left = (cell_types[cell_left] == 1) ? &data1 : &data2;
        ProblemData* data_right = (cell_types[cell_right] == 1) ? &data1 : &data2;
        if (data_left == data_right)
        {
            flux_calculators[i].reset(new InternalFluxCalculator(grid, *data_left, cell_left, cell_right));
        }
        else
        {
            flux_calculators[i].reset(
                new MergingFluxCalculator(grid, *data_left, *data_right, cell_left, cell_right, 1e-8));
        }
    }
    flux_calculators.back().reset(new OutflowFluxCalculator(grid, data1, grid.n_elements() - 1));
    AssemblerFlux assem(grid, data, flux_calculators);

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes(), data1.area0);

    AmgcMatrixSolver slv;
    NonstatGridSaver saver(grid1, generate_nodes_coo(gr1), "bebe");
    saver.new_time_step(0);
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(velocity, "velocity");

    size_t iter_max = 1;

    CsrMatrix block_u_transport;
    std::vector<double> coupling_flux_ua;
    std::vector<double> coupling_flux_u2;
    std::vector<double> coupling_flux_p;
    std::vector<double> tmp;
    double t = 0;
    // while (time < 0.25 - 1e-12){
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
        std::cout << "  P=" << p_inflow2(time) << std::endl;

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
        if (t > 0.0005 - 1e-6)
        {
            t = 0;
            saver.new_time_step(time);
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(velocity, "velocity");
        }

        std::vector<double> pressure = assem.pressure(area);
        monitor1.push_back(pressure[monitoring_node1]);
        monitor2.push_back(pressure[monitoring_node2]);
        monitor3.push_back(pressure[monitoring_node3]);
    }

    std::vector<double> pressure = assem.pressure(area);
    CHECK(pressure[20] == Approx(1521.4710994706).margin(1e-3));

    time_value_vtk(tau, monitor1, "monitor1.vtk");
    time_value_vtk(tau, monitor2, "monitor2.vtk");
    time_value_vtk(tau, monitor3, "monitor3.vtk");
}

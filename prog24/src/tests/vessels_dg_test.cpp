#include "bflow/graph_grid.hpp"
#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/vessel_graph.hpp"
#include "bflow/vtk.hpp"
#include "tests/nonlinsolver.hpp"
#include "bflow/fem_grid.hpp"
#include "bflow/flux_calculator.hpp"
#include "catch.hpp"
#include <fstream>

#include <iomanip>
using namespace bflow;

TEST_CASE("Single vessel, inviscid, explicit", "[single-vessel-inviscid-explicit]")
{
    ProblemData data;
    double time = 0;
    double L = 1.0;
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.01, 1);
    FemGrid grid(grid1);
    double tau = grid.h(0) / 100;
    std::vector<ElementBoundaryFluxes> upwind_fluxes(grid.n_elements());

    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    NonstatGridSaver saver(grid1, nodes_coo, "bububu");
    saver.new_time_step(0);

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes(), data.area0);
    std::vector<double> pressure(grid.n_nodes(), 0.0);
    saver.save_vtk_point_data(velocity, "velocity");
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(pressure, "pressure");

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
    auto qinput = [&time]()->double{
        return 1e-6 * exp(-1e4 * (time - 0.05) * (time - 0.05));
    };
    upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, qinput, 0));
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

    double t = 0.005;
    while (time < 0.2 - 1e-6)
    {
        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  Q=" << qinput() << std::endl;

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
            rhs_a[i] += tau * tran_a[i];
            rhs_u[i] += tau * tran_u[i];
        }
        // + coupling
        for (size_t ielem = 0; ielem < grid.n_elements(); ++ielem)
        {
            size_t node0 = grid.tab_elem_nodes(ielem)[0];
            size_t node1 = grid.tab_elem_nodes(ielem)[1];
            const auto& uf = upwind_fluxes[ielem];
            rhs_a[node0] += tau * data.flux_a(uf.upwind_area_x0, uf.upwind_velo_x0);
            rhs_a[node1] -= tau * data.flux_a(uf.upwind_area_x1, uf.upwind_velo_x1);
            rhs_u[node0] += tau * data.flux_u(uf.upwind_area_x0, uf.upwind_velo_x0);
            rhs_u[node1] -= tau * data.flux_u(uf.upwind_area_x1, uf.upwind_velo_x1);
        }
        // solve
        slv.solve(rhs_a, area);
        slv.solve(rhs_u, velocity);

        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            pressure[i] = data.pressure(area[i]);
        }
        t += tau;
        if (t - 1e-5< time )
        {
            saver.new_time_step(time);
            saver.save_vtk_point_data(velocity, "velocity");
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(pressure, "pressure");
            t += 0.005;
        }
    }
    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(24.3958).margin(1e-3));
}

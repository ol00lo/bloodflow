#include "bflow/matrix_solver.hpp"
#include "bflow/vessel_graph.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include "tests/test_utils.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace bflow;
namespace
{
std::mt19937 gen(0);
double ex(double x, double t = 0)
{
    if (t > 0)
        return ex(x - t);
    else
        return (x > 0 && x < 0.8) ? 1.0 : 0.0;
}
} // namespace

TEST_CASE("simple test vtk0", "[gridsaver-0]")
{
    // std::vector<std::vector<int>> node = {{0}, {0, 1}, {1, 2}, {2}};
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.2, 2);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");
    int i = string_count("first_try.vtk");
    CHECK(i == 42);
}

TEST_CASE("simple test vtk1", "[gridsaver-1]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2}};
    std::vector<double> ed = {1.4142135624, 1.0, 1.4142135624};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.8);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    int i = string_count("first_try.vtk");
    CHECK(i == 27);
}

TEST_CASE("simple test vtk2", "[gridsaver-2]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2, 3}, {3}};
    std::vector<double> ed = {1.4142135624, 1.0, 1.4142135624, 1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.5);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    std::vector<double> func1 = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
                                 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    std::vector<double> func2 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    std::vector<double> f1 = test_point_data(generate_points_coo(grid1, nodes_coo));
    vtk1.save_vtk_point_data(func1, "d1", "first_try.vtk");
    vtk1.save_vtk_cell_data(func2, "d2", "first_try.vtk");
    vtk1.save_vtk_point_data(func1, "d3", "first_try.vtk");
    vtk1.save_vtk_point_data(f1, "d4", "first_try.vtk");
    vtk1.save_vtk_cell_data(func2, "d5", "first_try.vtk");
    // vtk1.save_nonstat_vtkseries(5, 1, func2, "d1", "second_try.vtk");
    int i = string_count("first_try.vtk");
    CHECK(i == 139);
}

TEST_CASE("simple test vtk3", "[gridsaver-3]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2, 3}, {3}};
    std::vector<double> ed = {1.4142135624, 3.0, 2 * 1.4142135624, 2.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 1.0);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    int i = string_count("first_try.vtk");
    CHECK(i == 43);
}

TEST_CASE("simple test vtk4", "[gridsaver-4]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0, 4, 3}, {1, 5, 6, 7}, {2, 8, 9}, {3}, {4},
                                          {5},       {6},       {7},          {8},       {9}};
    std::vector<double> ed = {2.82, 4.0, 3.6, 3.16, 2.24, 3.6, 3.0, 2.24, 2.24, 3.6};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.2);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");
    std::vector<double> func1;
    std::vector<double> func2;
    for (int i = 0; i < grid1.n_nodes(); i++)
        func1.push_back(i);
    for (int i = 0; i < grid1.n_elem(); i++)
        func2.push_back(i);

    std::vector<double> f1 = test_point_data(generate_points_coo(grid1, nodes_coo));
    vtk1.save_vtk_point_data(func1, "d1", "first_try.vtk");
    vtk1.save_vtk_cell_data(func2, "d2", "first_try.vtk");
    vtk1.save_vtk_point_data(f1, "d3", "first_try.vtk");

    int i = string_count("first_try.vtk");
    CHECK(i == 1383);
}

TEST_CASE("nonstat", "[nonstat-1]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0, 4, 3}, {1, 5, 6, 7}, {2, 8, 9}, {3}, {4},
                                          {5},       {6},       {7},          {8},       {9}};
    std::vector<double> ed = {2.82, 4.0, 3.6, 3.16, 2.24, 3.6, 3.0, 2.24, 2.24, 3.6};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.2);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);

    double tend = 100;
    double tau = 1;
    NonstatGridSaver saver(grid1, nodes_coo, "second_try");
    std::vector<Point2> points_coo = generate_points_coo(grid1, nodes_coo);
    std::vector<double> for_cell;
    for (int i = 0; i < grid1.n_elem(); i++)
        for_cell.push_back(grid1.find_cell_length(i));

    for (double t = 0; t <= tend; t += tau)
    {
        std::vector<double> point_data = test_point_data(points_coo, t);
        std::vector<double> cell_data = test_cell_data(for_cell, t);
        saver.new_time_step(t);
        saver.save_vtk_point_data(point_data, "d1");
        saver.save_vtk_cell_data(cell_data, "d2");
    }
    int i = string_count("second_try.vtk.series");
    CHECK(i == 106);
}

TEST_CASE("test dgrid", "[dgrid1]")
{
    std::vector<std::vector<int>> node = {{0, 1}, {0, 2}, {1, 3, 4}, {2}, {3}, {4}};
    std::vector<double> ed = {1.0, 1.0, 1.0, 1.0, 1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.5, 2);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    NonstatGridSaver saver(grid1, nodes_coo, "second_try11");
    for (double t = 0; t <= 10; t += 0.2)
    {
        std::vector<double> point_data;
        for (int i = 0; i < grid1.n_nodes(); i++)
            point_data.push_back(std::sin(i * t));
        saver.new_time_step(t);
        saver.save_vtk_point_data(point_data, "d1");
    }
    int i = string_count("second_try11.vtk.series");
    CHECK(i == 56);
}

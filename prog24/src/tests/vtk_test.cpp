#include "bflow/vessel_graph.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include <fstream>
#include <iostream>
#include <vector>

using namespace bflow;

TEST_CASE("simple test vtk0", "[gridsaver-0]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.25);
    std::vector<Point2> nodes_coo = grid1.generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    char *str = new char[1024];
    int i = 0;
    std::ifstream base("first_try.vtk");
    while (!base.eof())
    {
        base.getline(str, 1024, '\n');
        i++;
    }
    CHECK(i == 21);
}

TEST_CASE("simple test vtk1", "[gridsaver-1]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2}};
    std::vector<double> ed = {1.4142135624, 1.0, 1.4142135624};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.8);
    std::vector<Point2> nodes_coo = grid1.generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    char *str = new char[1024];
    int i = 0;
    std::ifstream base("first_try.vtk");
    while (!base.eof())
    {
        base.getline(str, 1024, '\n');
        i++;
    }
    CHECK(i == 24);
}

TEST_CASE("simple test vtk2", "[gridsaver-2]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2, 3}, {3}};
    std::vector<double> ed = {1.4142135624, 1.0, 1.4142135624, 1.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 1.0);
    std::vector<Point2> nodes_coo = grid1.generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    char *str = new char[1024];
    int i = 0;
    std::ifstream base("first_try.vtk");
    while (!base.eof())
    {
        base.getline(str, 1024, '\n');
        i++;
    }
    CHECK(i == 21);
}

TEST_CASE("simple test vtk3", "[gridsaver-3]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2, 3}, {3}};
    std::vector<double> ed = {1.4142135624, 3.0, 2 * 1.4142135624, 2.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 1.0);
    std::vector<Point2> nodes_coo = grid1.generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    char *str = new char[1024];
    int i = 0;
    std::ifstream base("first_try.vtk");
    while (!base.eof())
    {
        base.getline(str, 1024, '\n');
        i++;
    }
    CHECK(i == 36);
}

TEST_CASE("simple test vtk4", "[gridsaver-4]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0, 4, 3}, {1, 5, 6, 7}, {2, 8, 9}, {3}, {4},
                                          {5},       {6},       {7},          {8},       {9}};
    std::vector<double> ed = {2.82, 4.0, 3.6, 3.16, 2.24, 3.6, 3.0, 2.24, 2.24, 3.6};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.2);
    std::vector<Point2> nodes_coo = grid1.generate_nodes_coo(gr1);
    GridSaver vtk1(grid1, nodes_coo);
    vtk1.save_area("first_try.vtk");

    char *str = new char[1024];
    int i = 0;
    std::ifstream base("first_try.vtk");
    while (!base.eof())
    {
        base.getline(str, 1024, '\n');
        i++;
    }
    CHECK(i == 465);
}
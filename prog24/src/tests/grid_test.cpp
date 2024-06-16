#include "bflow/graph_grid.hpp"
#include "bflow/vessel_graph.hpp"
#include "catch.hpp"
#include <map>
#include <vector>

using namespace bflow;

TEST_CASE("simple test grid", "[simple_grid_test]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {0.5};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.61);
    CHECK(grid1.n_points() == 2);
    // CHECK(grid1.n_cells() == 1);
    CHECK(grid1.points_by_edge(0)[1] == 1);
    // CHECK(grid1.find_cell_length(0) == 0.5);
    node = {{0}, {0}};
    ed = {3.0};
    VesselGraph gr2(node, ed);
    GraphGrid grid2(gr2, 0.5);
    CHECK(grid2.n_points() == 7);
    // CHECK(grid2.n_cells() == 6);
    CHECK(grid2.points_by_edge(0)[1] == 2);
    // CHECK(grid2.find_cell_length(0) == 0.5);
    node = {{0}, {0, 1, 2}, {1}, {2}};
    ed = {3.0, 0.5, 10.0};
    VesselGraph gr3(node, ed);
    GraphGrid grid3(gr3, 0.6);
    CHECK(grid3.n_points() == 24);
    // CHECK(grid3.n_cells() == 23);
    CHECK(grid3.points_by_edge(0)[1] == 4);
    // CHECK(grid3.find_cell_length(0) == 0.6);
}

TEST_CASE("test grid", "[grid_test]")
{
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0, 4, 3}, {1, 5, 6, 7}, {2, 8, 9}, {3}, {4},
                                          {5},       {6},       {7},          {8},       {9}};
    std::vector<double> ed = {2.82, 4.0, 3.6, 3.16, 2.24, 3.6, 3.0, 2.24, 2.24, 3.6};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.61);
    CHECK(grid1.tab_cell_point(2)[0] == 12);
    CHECK(grid1.tab_cell_point(2)[1] == 13);
    CHECK(grid1.tab_point_cell(0).size() == 3);
    CHECK(grid1.tab_point_cell(2).size() == 4);
    CHECK(grid1.tab_point_cell(0)[1] == 5);
    CHECK(grid1.tab_point_cell(12).size() == 2);
    // CHECK(grid1.find_edge_by_cell(7) == 1);
    // CHECK_THROWS(grid1.find_edge_by_cell(55));
}

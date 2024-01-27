#include "bflow/vessel_graph.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include <vector>

using namespace bflow;

TEST_CASE("simple test vtk", "[simple_vtk_test]")
{
    // std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2}};
    // std::vector<double> ed = {1.0, 1.0, 1.0};
    // std::vector<std::vector<int>> node = {{0, 1, 2}, {0}, {1}, {2, 3, 4}, {3}, {4}};
    // std::vector<double> ed = {1.0, 1.0, 1.0, 0.5, 0.5};
    std::vector<std::vector<int>> node = {{0, 1, 2}, {0, 4, 3}, {1, 5, 6, 7}, {2, 8, 9}, {3}, {4},
                                          {5},       {6},       {7},          {8},       {9}};
    std::vector<double> ed = {2.82, 4.0, 3.6, 3.16, 2.24, 3.6, 3.0, 2.24, 2.24, 3.6};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.5);
    GridSaver vtk1(grid1);
    vtk1.points_save("first_try.vtk");
}
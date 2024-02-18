#ifndef GRAPH_GRID_HPP
#define GRAPH_GRID_HPP
#include "vessel_graph.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace bflow
{
struct Point2
{
    double x;
    double y;
};

class GraphGrid
{
public:
    GraphGrid(const VesselGraph& graph, double h);
    int n_points() const;
    int n_cells() const;
    int n_edges() const;
    std::vector<int> tab_point_cell(int point) const;
    std::array<int, 2> tab_cell_point(int cell) const;
    int find_edge_by_cell(int cell) const;
    std::vector<int> points_by_edge(int edge) const;
    double find_cell_length(int cell) const;
    std::array<int, 2> find_node_by_edge(int edge) const;

private:
    std::vector<std::vector<int>> _points;
    std::vector<double> _cells;

    std::vector<std::vector<int>> _point_cells;
    std::vector<std::array<int, 2>> _cell_points;
    std::vector<int> _cell_edges;
    int _n_points;
    std::vector<std::array<int, 2>> _edge_nodes;
    // int _n_cells;
};
} // namespace bflow
#endif
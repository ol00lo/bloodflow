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
    GraphGrid(const VesselGraph& graph, double h, int nadd=1);
    int n_points() const;
    int n_nodes() const
    {
        return _n_nodes;
    }
    int n_cells() const;
    int n_edges() const;
    std::vector<int> tab_point_cell(int point) const;
    std::array<int, 2> tab_cell_point(int cell) const;
    int find_edge_by_cell(int cell) const;
    std::vector<int> points_by_edge(int edge) const;
    std::vector<int> nodes_by_edge(int edge) const
    {
        return _nodes_by_edge[edge];
    }
    double find_cell_length(int cell) const;
    std::array<int, 2> find_node_by_edge(int edge) const;
    const int _nadd;
    std::array<int, 2> node_by_cell(int cell) const
    {
        return _cells[cell];
    };

private:
    std::vector<std::vector<int>> _points;
    std::vector<double> _cellslen;
    std::vector<std::array<int, 2>> _cells;
    std::vector<std::vector<int>> _point_cells;
    std::vector<std::array<int, 2>> _cell_points;
    std::vector<int> _cell_edges;
    int _n_points;
    int _n_nodes;
    std::vector<std::array<int, 2>> _edge_points;
    std::vector<std::vector<int>> _nodes_by_edge;
};
} // namespace bflow
#endif
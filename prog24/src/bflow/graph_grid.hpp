#ifndef GRAPH_GRID_HPP
#define GRAPH_GRID_HPP
#include "vessel_graph.hpp"
#include "bflow/matrix.hpp"
#include "bflow/macros.hpp"
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
    int n_nodes() const;
    int n_elem() const;
    int n_edges() const;
    std::vector<int> tab_point_cell(int point) const;
    std::array<int, 2> tab_cell_point(int cell) const;
    int find_edge_by_cell(int cell) const;
    std::vector<int> points_by_edge(int edge) const;
    std::vector<int> nodes_by_edge(int edge) const;
    double find_cell_length(int cell) const;
    std::array<int, 2> find_node_by_edge(int edge) const;
    std::array<int, 2> find_point_by_edge(int edge) const;
    const int _power;
    std::array<int, 2> node_by_cell(int cell) const;
    std::vector<std::array<int, 2>> cells() const;
    std::vector<int> tab_point_nodes(int ipoint) const
    {
        return (*_point_nodes.find(ipoint)).second;
    }

private:
    std::vector<std::vector<int>> _points_by_edge;
    std::vector<double> _cellslen;
    std::vector<std::array<int, 2>> _cells;
    std::vector<std::vector<int>> _point_cells;
    std::vector<std::array<int, 2>> _cell_points;
    std::vector<int> _cell_edges;
    int _n_points;
    int _n_nodes;
    std::vector<int> _local_nodes;
    std::vector<std::array<int, 2>> _edge_points;
    std::vector<std::vector<int>> _nodes_by_edge;
    std::vector<std::array<int, 2>> _bound_points;
    std::map<int, std::vector<int>> _point_nodes;

};
} // namespace bflow
#endif
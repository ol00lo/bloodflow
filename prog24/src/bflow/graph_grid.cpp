#include "graph_grid.hpp"
#define PI 3.1415926

using namespace bflow;

GraphGrid::GraphGrid(const VesselGraph& graph, double h)
{
    _points.resize(graph.n_edges());
    int n_nodes = graph.n_nodes();

    for (int edge = 0; edge < graph.n_edges(); ++edge)
    {

        std::array<int, 2> boundary = graph.tab_edge_node(edge);
        _edge_nodes.push_back(boundary);
        double len = graph.find_length(edge);
        int m_cells = std::round(len / h);
        if (len < 0.5 * h)
        {
            m_cells = 1;
        }
        double h_cell = 1.0 * len / m_cells;

        _points[edge].push_back(boundary[0]);
        for (size_t cells = 1; cells < m_cells; ++cells)
        {
            _points[edge].push_back(n_nodes);
            _cells.push_back(h_cell * cells);
            n_nodes++;
        }
        _points[edge].push_back(boundary[1]);
        _cells.push_back(len);
    }

    std::set<int> m_points;
    for (size_t i = 0; i < _points.size(); ++i)
    for (size_t j = 0; j < _points[i].size(); ++j)
    {
        m_points.insert(_points[i][j]);
    }
    _n_points = m_points.size();

    for (size_t i = 0; i < _points.size(); ++i)
    for (size_t j = 1; j < _points[i].size(); ++j)
    {
        std::array<int, 2> points;
        points[0] = _points[i][j - 1];
        points[1] = _points[i][j];
        _cell_points.push_back(points);
    }

    _point_cells.resize(_n_points);
    for (int cell = 0; cell < _cell_points.size(); ++cell)
    {
        _point_cells[_cell_points[cell][0]].push_back(cell);
        _point_cells[_cell_points[cell][1]].push_back(cell);
    }

    for (size_t i = 0; i < _points.size(); ++i)
    for (size_t j = 1; j < _points[i].size(); ++j)
    {
        _cell_edges.push_back(i);
    }
}

int GraphGrid::n_points() const
{
    return _n_points;
}

int GraphGrid::n_cells() const
{
    return _cells.size();
}

std::array<int, 2> GraphGrid::tab_cell_point(int cell) const
{
    return _cell_points.at(cell);
}

std::vector<int> GraphGrid::tab_point_cell(int point) const
{
    return _point_cells.at(point);
}

int GraphGrid::find_edge_by_cell(int cell) const
{
    return _cell_edges.at(cell);
}

std::vector<int> GraphGrid::points_by_edge(int edge) const
{
    return _points.at(edge);
}

double GraphGrid::find_cell_length(int cell) const
{
    return _cells.at(cell);
}

int GraphGrid::n_edges() const
{
    return _edge_nodes.size();
}

std::array<int, 2> GraphGrid::find_node_by_edge(int edge) const
{
    return _edge_nodes.at(edge);
}
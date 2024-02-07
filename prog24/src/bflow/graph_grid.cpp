#include "graph_grid.hpp"
#include <algorithm>
#define PI 3.1415926

using namespace bflow;

GraphGrid::GraphGrid(const VesselGraph &graph, double h)
{
    _points.resize(graph.n_edges());
    int n_nodes = graph.n_nodes();

    for (int edge = 0; edge < graph.n_edges(); ++edge)
    {

        std::array<int, 2> boundary = graph.tab_edge_node(edge);
        double len = graph.find_length(edge);
        int n_cells = std::round(len / h);
        if (len < 0.5 * h)
        {
            n_cells = 1;
        }
        double h_cell = 1.0 * len / n_cells;

        _points[edge].push_back(boundary[0]);
        for (size_t cells = 1; cells < n_cells; ++cells)
        {
            _points[edge].push_back(n_nodes);
            _cells.push_back(h_cell * cells);
            n_nodes++;
        }
        _points[edge].push_back(boundary[1]);
        _cells.push_back(len);
    }

    std::set<int> n_points;
    for (size_t i = 0; i < _points.size(); ++i)
        for (size_t j = 0; j < _points[i].size(); ++j)
        {
            n_points.insert(_points[i][j]);
        }
    _n_points = n_points.size();

    //_n_cells = 0;
    // for (size_t i = 0; i < _cells.size(); ++i) {
    //	_n_cells += _cells[i].size();
    //}

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

    int i_cell = 0;
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
    if (cell >= _cells.size())
    {
        throw std::runtime_error("Cell out of grid");
    }
    return _cell_points[cell];
}

std::vector<int> GraphGrid::tab_point_cell(int point) const
{
    if (point > _n_points)
    {
        throw std::runtime_error("Point out of grid");
    }
    return _point_cells[point];
}

int GraphGrid::find_edge_by_cell(int cell) const
{
    if (cell >= _cells.size())
    {
        throw std::runtime_error("Cell out of grid");
    }
    return _cell_edges[cell];
}

std::vector<int> GraphGrid::points_by_edge(int edge) const
{
    if (edge > _points.size())
    {
        throw std::runtime_error("Cell out of grid");
    }
    return _points[edge];
}

double GraphGrid::find_cell_length(int cell) const
{
    if (cell >= _cells.size())
    {
        throw std::runtime_error("Cell out of grid");
    }
    return _cells[cell];
}

std::vector<Point2> GraphGrid::generate_nodes_coo(const VesselGraph &graph) const
{
    std::vector<Point2> nodes;
    Point2 old_point = {0.0, 0.0};
    nodes.push_back(old_point);

    int edges = -1;
    for (int node = 0; node < graph.n_nodes(); ++node)
    {
        old_point = nodes[node];
        std::vector<int> nedges = graph.tab_node_edge(node);
        std::sort(nedges.begin(), nedges.end());
        int size = nedges.size();

        int i = 1;
        for (int edge = 0; edge < size; ++edge)
        {
            if (nedges[edge] > edges)
            {
                edges++;
                Point2 new_point;
                double len = graph.find_length(edges);
                double angle = PI - PI * i / size;
                if (node == 0)
                    angle = PI - PI * i / (size + 1);
                new_point.x = old_point.x + len * std::cos(angle);
                new_point.y = old_point.y + len * std::sin(angle);
                nodes.push_back(new_point);
                i++;
            }
        }
    }

    return nodes;
}
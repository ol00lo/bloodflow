#include "graph_grid.hpp"
#define PI 3.1415926

using namespace bflow;

GraphGrid::GraphGrid(const VesselGraph& graph, double h, int nadd) : n_midnodes(nadd)
{
    _points_by_edge.resize(graph.n_edges());
    _nodes_by_edge.resize(graph.n_edges());
    int n_points = graph.n_nodes();
    int node = 0;

    for (int edge = 0; edge < graph.n_edges(); ++edge)
    {
        std::array<int, 2> boundary = graph.tab_edge_node(edge);
        _bound_points.push_back(boundary);
        double len = graph.find_length(edge);
        int m_cells = std::round(len / h);
        if (len < 0.5 * h)
        {
            m_cells = 1;
        }
        double h_cell = 1.0 * len / m_cells;

        _nodes_by_edge[edge].push_back(node);
        _points_by_edge[edge].push_back(boundary[0]);
        _point_nodes.insert({boundary[0], {node}});
        node++;
        for (size_t cells = 1; cells < m_cells; ++cells)
        {
            _points_by_edge[edge].push_back(n_points);
            _nodes_by_edge[edge].push_back(node);
            node++;
            _nodes_by_edge[edge].push_back(node);
            _point_nodes.insert({n_points, {node-1, node}});
            node++;
            //_cellslen.push_back(h_cell * cells);
            _cellslen.push_back(h_cell);
            n_points++;
        }
        _points_by_edge[edge].push_back(boundary[1]);
        _nodes_by_edge[edge].push_back(node);
        _point_nodes.insert({boundary[1], {node}});
        node++;
        _cellslen.push_back(h_cell);

        _edge_points.push_back({_nodes_by_edge[edge][0], _nodes_by_edge[edge].back()});
    }
    _n_nodes = node;

    for (int edge = 0; edge < _nodes_by_edge.size(); edge++)
    {
        int size = _nodes_by_edge[edge].size();
        for (int nod = 0; nod < size; nod += 2)
        {
            int before = _nodes_by_edge[edge][nod];
            int after = _nodes_by_edge[edge][nod + 1];
            if (n_midnodes > 1)
            {
                for (int i = 0; i < n_midnodes - 1; i++)
                {
                    _nodes_by_edge[edge].push_back(node);
                    after = node;
                    _cells.push_back({before, after});
                    before = node;
                    node++;
                }
            }
            after = _nodes_by_edge[edge][nod + 1];
            _cells.push_back({before, after});
        }
    }

    _n_nodes = node;

    std::set<int> m_points;
    for (size_t i = 0; i < _points_by_edge.size(); ++i)
        for (size_t j = 0; j < _points_by_edge[i].size(); ++j)
        {
            m_points.insert(_points_by_edge[i][j]);
        }
    _n_points = m_points.size();

    for (size_t i = 0; i < _points_by_edge.size(); ++i)
        for (size_t j = 1; j < _points_by_edge[i].size(); ++j)
        {
            std::array<int, 2> points;
            points[0] = _points_by_edge[i][j - 1];
            points[1] = _points_by_edge[i][j];
            _cell_points.push_back(points);
        }

    _point_cells.resize(_n_points);
    for (int cell = 0; cell < _cell_points.size(); ++cell)
    {
        _point_cells[_cell_points[cell][0]].push_back(cell);
        _point_cells[_cell_points[cell][1]].push_back(cell);
    }

}

int GraphGrid::n_points() const
{
    return _n_points;
}

 int GraphGrid::n_nodes() const
{
    return _n_nodes;
}

int GraphGrid::n_elem() const
{
    return _cellslen.size();
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
    // return 0;
}

std::vector<int> GraphGrid::points_by_edge(int edge) const
{
    return _points_by_edge.at(edge);
}

std::vector<int> GraphGrid::nodes_by_edge(int edge) const
{
    return _nodes_by_edge[edge];
}

double GraphGrid::find_cell_length(int cell) const
{
    return _cellslen.at(cell);
}

int GraphGrid::n_edges() const
{
    return _edge_points.size();
}

std::array<int, 2> GraphGrid::find_node_by_edge(int edge) const
{
    return _edge_points.at(edge);
}

std::array<int, 2> GraphGrid::find_points_by_edge(int edge) const
{
    return _bound_points[edge];
}

std::array<int, 2> GraphGrid::node_by_cell(int cell) const
{
    return _cells[cell];
};

 std::vector<std::array<int, 2>> GraphGrid::cells() const
{
    return _cells;
}
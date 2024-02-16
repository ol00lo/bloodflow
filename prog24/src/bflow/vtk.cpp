#include "vtk.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#define PI 3.1415926

using namespace bflow;

std::vector<Point2> bflow::generate_nodes_coo(const VesselGraph& graph)
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
                double angle = (node == 0) ? (PI - PI * i / (size + 1)) : (PI - PI * i / size);
                new_point.x = old_point.x + len * std::cos(angle);
                new_point.y = old_point.y + len * std::sin(angle);
                nodes.push_back(new_point);
                i++;
            }
        }
    }

    return nodes;
}

GridSaver::GridSaver(const GraphGrid& grid, const std::vector<Point2>& nodes_coo)
{
    for (int iedge = 0; iedge < grid.n_edges(); iedge++)
    {
        std::vector<int> points_by_edge = grid.points_by_edge(iedge);
        std::array<int, 2> edge_nodes = grid.find_node_by_edge(iedge);
        _points.push_back(nodes_coo[edge_nodes[0]]);
        int n_edge_cells = points_by_edge.size() - 1;
        int i = 1;
        for (int ipoint = 0; ipoint < n_edge_cells - 1; ipoint++)
        {
            Point2 new_point;
            double w = (double)i / n_edge_cells;
            new_point.x = (1 - w) * nodes_coo[edge_nodes[0]].x + w * nodes_coo[edge_nodes[1]].x;
            new_point.y = (1 - w) * nodes_coo[edge_nodes[0]].y + w * nodes_coo[edge_nodes[1]].y;
            _points.push_back(new_point);
            i++;
            int a = _points.size();
            _cells.push_back({a - 2, a - 1});
        }
        _points.push_back(nodes_coo[edge_nodes[1]]);
        int a = _points.size();
        _cells.push_back({a - 2, a- 1});
    }
}

void GridSaver::save_area(std::string filename) const
{
    std::ofstream fs(filename);
    fs << "# vtk DataFile Version 3.0" << std::endl;
    fs << "Unstructured grid" << std::endl;
    fs << "ASCII" << std::endl;
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;

    fs << "POINTS " << _points.size() << " double" << std::endl;
    for (int ipoint = 0; ipoint < _points.size(); ++ipoint)
    {
        Point2 point = _points[ipoint];
        fs << point.x << " " << point.y << " " << 0.0 << std::endl;
    }

    fs << "CELLS  " << _cells.size() << "   " << 3 * _cells.size() << std::endl;
    for (size_t i = 0; i < _cells.size(); ++i)
    {
        fs << 2 << " " << _cells[i][0] << " " << _cells[i][1] << std::endl;
    }

    fs << "CELL_TYPES  " << _cells.size() << std::endl;
    for (size_t i = 0; i < _cells.size(); ++i)
        fs << 3 << std::endl;
}
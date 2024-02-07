#include "vtk.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#define PI 3.1415926

using namespace bflow;
//
// GridSaver::GridSaver(const GraphGrid &grid)
//{
//    _points.resize(grid.n_points());
//
//    Point2 zero;
//    zero.x = 0.0;
//    zero.y = 0.0;
//    _points[0] = zero;
//
//    double angle = 0.0;
//    double len = 0.0;
//    int edge = 0;
//
//    for (int cell = 0; cell < grid.n_cells(); ++cell)
//    {
//        std::array<int, 2> ipoints = grid.tab_cell_point(cell);
//        Point2 c_point_old = _points[ipoints[0]];
//        c_point_old.y = _points[grid.tab_cell_point(cell)[0]].y;
//
//        std::vector<int> for_cells(ipoints.begin(), ipoints.end());
//        _cells.push_back(for_cells);
//
//        if (grid.tab_cell_point(cell)[0] == 0)
//        {
//            edge++;
//            len = grid.find_cell_length(cell);
//            angle = PI / grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size();
//        }
//        else if (grid.find_edge_by_cell(cell) != grid.find_edge_by_cell(cell - 1))
//        {
//            edge++;
//            len = grid.find_cell_length(cell);
//            angle = PI / grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size();
//            if (edge >= grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size())
//            {
//                edge = 1;
//            }
//        }
//
//        if (grid.tab_cell_point(cell)[0] == 0)
//        {
//            angle = PI / (grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size() + 1);
//        }
//
//        Point2 c_point_new;
//        c_point_new.x = c_point_old.x + len * std::cos(PI - edge * angle);
//        c_point_new.y = c_point_old.y + len * std::sin(PI - edge * angle);
//
//        int point = grid.tab_cell_point(cell)[1];
//        _points[point] = c_point_new;
//    }
//}

GridSaver::GridSaver(const GraphGrid &grid, const std::vector<Point2> &nodes_coo)
{
    double angle = 0.0;
    int edge = 0;
    int cell = 0;
    _points.push_back(nodes_coo[0]);
    for (int node = 0; node < nodes_coo.size(); ++node)
    {
        int size = grid.tab_point_cell(node).size() - 1;
        if (node == 0)
            size++;
        if (size > 0)
        {
            int for_cells0 = 0;
            for (auto it = _points.begin(); it < _points.end(); it++)
            {
                Point2 p = *it;
                if (p.x == nodes_coo[node].x && p.y == nodes_coo[node].y)
                {
                    for_cells0 = it - _points.begin();
                }
            }
            int i = 1;
            for (int ed = 0; ed < size; ++ed)
            {
                angle = PI - PI * i / (size + 1);
                Point2 cpoint_old = nodes_coo[node];

                Point2 cpoint_new;

                std::vector<int> points_by_edge = grid.points_by_edge(edge);

                std::array<int, 2> for_cells;
                for_cells[0] = for_cells0;
                for_cells[1] = _points.size();

                _cells.push_back(for_cells);

                for (int j = 0; j < points_by_edge.size() - 2; ++j)
                {
                    std::array<int, 2> for_cells;
                    for_cells[0] = _points.size();
                    for_cells[1] = _points.size() + 1;
                    _cells.push_back(for_cells);
                    double len = grid.find_cell_length(cell);
                    cpoint_new.x = cpoint_old.x + len * std::cos(angle);
                    cpoint_new.y = cpoint_old.y + len * std::sin(angle);
                    _points.push_back(cpoint_new);
                    cell++;
                }

                cell++;
                edge++;
                i++;
                _points.push_back(nodes_coo[edge]);
            }
        }
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
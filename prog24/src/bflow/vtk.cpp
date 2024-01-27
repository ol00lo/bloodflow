#include "vtk.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#define PI 3.1415926

using namespace bflow;

GridSaver::GridSaver(const GraphGrid &grid)
{
    _points.resize(grid.n_points());

    Point zero;
    zero.x = 0.0;
    zero.y = 0.0;
    _points[0] = zero;

    double angle = 0.0;
    double len = 0.0;
    int edge = 0;

    for (int cell = 0; cell < grid.n_cells(); ++cell)
    {
        Point c_point_old;
        c_point_old.x = _points[grid.tab_cell_point(cell)[0]].x;
        c_point_old.y = _points[grid.tab_cell_point(cell)[0]].y;

        std::vector<int> for_cells;
        for_cells.push_back(grid.tab_cell_point(cell)[0]);
        for_cells.push_back(grid.tab_cell_point(cell)[1]);
        _cells.push_back(for_cells);

        if (grid.tab_cell_point(cell)[0] == 0)
        {
            edge++;
            len = grid.find_cell_length(cell);
            angle = PI / grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size();
        }
        else if (grid.find_edge_by_cell(cell) != grid.find_edge_by_cell(cell - 1))
        {
            edge++;
            len = grid.find_cell_length(cell);
            angle = PI / grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size();
            if (edge >= grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size())
            {
                edge = 1;
            }
        }

        if (grid.tab_cell_point(cell)[0] == 0)
        {
            angle = PI / (grid.tab_point_cell(grid.tab_cell_point(cell)[0]).size() + 1);
        }

        Point c_point_new;
        c_point_new.x = c_point_old.x + len * std::cos(PI - edge * angle);
        c_point_new.y = c_point_old.y + len * std::sin(PI - edge * angle);

        int point = grid.tab_cell_point(cell)[1];
        _points[point] = c_point_new;
    }
}

void GridSaver::points_save(std::string filename) const
{
    std::ofstream fs(filename);
    fs << "# vtk DataFile Version 3.0" << std::endl;
    fs << "Unstructured grid" << std::endl;
    fs << "ASCII" << std::endl;
    fs << "DATASET UNSTRUCTURED_GRID" << std::endl;

    fs << "POINTS " << _points.size() << " double" << std::endl;
    for (int ipoint = 0; ipoint < _points.size(); ++ipoint)
    {
        Point point = _points[ipoint];
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
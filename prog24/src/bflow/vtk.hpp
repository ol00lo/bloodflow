#ifndef VTK_HPP
#define VTK_HPP

#include "graph_grid.hpp"

struct Point
{
    double x;
    double y;
};

namespace bflow
{
class GridSaver
{
public:
    GridSaver(const GraphGrid &grid);
    void points_save(std::string filename) const;

private:
    std::vector<Point> _points;
    std::vector<std::vector<int>> _cells;
};
} // namespace bflow

#endif
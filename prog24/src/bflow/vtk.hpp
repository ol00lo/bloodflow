#ifndef VTK_HPP
#define VTK_HPP

#include "graph_grid.hpp"

namespace bflow
{
std::vector<Point2> generate_nodes_coo(const VesselGraph& graph);
class GridSaver
{
public:
    GridSaver(const GraphGrid& grid, const std::vector<Point2>& nodes_coo);
    void save_area(std::string filename) const;

private:
    std::vector<Point2> _points;
    std::vector<std::array<int, 2>> _cells;
};
} // namespace bflow

#endif
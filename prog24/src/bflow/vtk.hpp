#ifndef VTK_HPP
#define VTK_HPP

#include "graph_grid.hpp"

namespace bflow
{
std::vector<Point2> generate_nodes_coo(const VesselGraph& graph);
std::vector<Point2> generate_points_coo(const GraphGrid& grid, const std::vector<Point2>& nodes_coo);

class GridSaver
{
public:
    GridSaver(const GraphGrid& grid, const std::vector<Point2>& nodes_coo);
    void save_area(std::string filename) const;
    void save_vtk_point_data(const std::vector<double>& vertex_data, std::string dataname, std::string filename) const;
    void save_vtk_cell_data(const std::vector<double>& cell_data, std::string dataname, std::string filename) const;

private:
    std::vector<Point2> _points;
    std::vector<std::array<int, 2>> _cells;
};


class NonstatGridSaver
{
public:
    NonstatGridSaver(const GraphGrid& grid, const std::vector<Point2>& nodes_coo, std::string filename);
    void new_time_step(double t);
    void save_vtk_point_data(const std::vector<double>& data, std::string data_name) const;
    void save_vtk_cell_data(const std::vector<double>& data, std::string data_name) const;
private:
    GridSaver _vtk;
    std::vector<double> _times;
    std::string _file_name;
    std::string _series_name;
    std::vector<std::string> _files;
    int _files_count = 0;
    void add_in_series() const;
    std::string print_files() const;
};


} // namespace bflow

#endif
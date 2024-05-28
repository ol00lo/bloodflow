#include "vtk.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#define PI 3.1415926

using namespace bflow;

namespace
{
std::vector<Point2> find_coo(Point2 x, Point2 y, int col)
{
    std::vector<Point2> res;
    int i = 1;
    for (int ipoint = 0; ipoint < col - 1; ipoint++)
    {
        Point2 new_point;
        double w = (double)i / col;
        new_point.x = (1 - w) * x.x + w * y.x;
        new_point.y = (1 - w) * x.y + w * y.y;
        res.push_back(new_point);
        i++;
    }
    return res;
}

std::vector<std::string> to_strings(std::string file_name)
{
    std::vector<std::string> res;
    std::ifstream in(file_name);
    std::string str;
    while (std::getline(in, str))
    {
        res.push_back(str);
    }
    return res;
}

void write_data(std::fstream& outfile, std::string dataname, int size, const std::vector<double>& data)
{
    outfile << "SCALARS " << dataname << " double 1" << std::endl;
    outfile << "LOOKUP_TABLE default" << std::endl;

    for (size_t i = 0; i < size; ++i)
        outfile << data[i] << std::endl;
}

void write_new_data(std::string filename, std::string type, std::string dataname, int size,
                    const std::vector<double>& data)
{
    std::fstream fs(filename, std::ios::app);
    fs << type << " " << size << std::endl;
    write_data(fs, dataname, size, data);
}

void write_not_new_data(std::string filename, int lines_count, std::string dataname, int size,
                        const std::vector<double>& data)
{
    std::vector<std::string> help;
    help = to_strings(filename);
    if (help.size() <= lines_count + size + 2)
    {
        std::fstream fs(filename, std::ios::app);
        write_data(fs, dataname, size, data);
    }
    else
    {
        std::fstream fs(filename);
        for (int line = 0; line < help.size(); line++)
        {
            if (line == lines_count + size + 2)
            {
                write_data(fs, dataname, size, data);
            }
            fs << help[line] << std::endl;
        }
    }
}

int find_data_in_file(std::string filename, std::string data)
{
    int size = data.size();
    std::ifstream fse(filename);
    std::string str;
    bool found = false;
    int sum_str = 0;
    while (std::getline(fse, str))
    {
        sum_str++;
        if (str.substr(0, size) == data)
        {
            found = true;
            break;
        }
    }
    if (found == true)
        return sum_str;
    else
        return -1;
}
} // namespace

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

std::vector<Point2> bflow::generate_points_coo(const GraphGrid& grid, const std::vector<Point2>& nodes_coo)
{
    std::vector<Point2> points_coo(grid.n_nodes());
    for (int iedge = 0; iedge < grid.n_edges(); iedge++)
    {
        std::vector<Point2> added;
        std::vector<int> nodes_by_edge = grid.nodes_by_edge(iedge);
        std::vector<int> points_by_edge = grid.points_by_edge(iedge);
        std::array<int, 2> edge_nodes = grid.find_node_by_edge(iedge);
        std::array<int, 2> p_nodes = grid.find_point_by_edge(iedge);
        points_coo[edge_nodes[0]]=nodes_coo[p_nodes[0]];
        int n_edge_cells = points_by_edge.size() - 1;
        int i = 1;
        auto ipoint = nodes_by_edge.begin()+1;
        while (*ipoint!=edge_nodes.back())
        {
            Point2 new_point;
            double w = (double)i / n_edge_cells;
            new_point.x = (1 - w) * nodes_coo[p_nodes[0]].x + w * nodes_coo[p_nodes[1]].x;
            new_point.y = (1 - w) * nodes_coo[p_nodes[0]].y + w * nodes_coo[p_nodes[1]].y;
            points_coo[*ipoint]=new_point;
            points_coo[*(ipoint+1)]=new_point;
            auto x =
                find_coo(points_coo[*(ipoint - 1)], points_coo[*ipoint], grid._power);
            for (int j = 0; j < grid._power-1; j++)
            {
                added.push_back(x[j]);
            }
            ipoint+=2;
            i++;
        }
        points_coo[edge_nodes[1]] = nodes_coo[p_nodes[1]];
        auto x = find_coo(points_coo[*(ipoint - 1)], points_coo[*ipoint], grid._power);
        for (int j = 0; j < grid._power - 1; j++)
        {
            added.push_back(x[j]);
        }
        for (auto x : added)
        {
            ipoint++;
            points_coo[*ipoint]=x;
        }
    }
    return points_coo;
}

GridSaver::GridSaver(const GraphGrid& grid, const std::vector<Point2>& nodes_coo) : _cells(grid.cells())
{
    _points = generate_points_coo(grid, nodes_coo);
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

void GridSaver::save_vtk_point_data(const std::vector<double>& vertex_data, std::string dataname,
                                    std::string filename) const
{
    if (vertex_data.size() != _points.size())
        throw std::runtime_error("incorrect number of vertex data");

    int str_to_data = find_data_in_file(filename, "POINT_DATA");
    if (str_to_data == -1)
        write_new_data(filename, "POINT_DATA", dataname, _points.size(), vertex_data);
    else
        write_not_new_data(filename, str_to_data, dataname, _points.size(), vertex_data);
}

void GridSaver::save_vtk_cell_data(const std::vector<double>& cell_data, std::string dataname,
                                   std::string filename) const
{
    if (cell_data.size() != _cells.size())
        throw std::runtime_error("incorrect number of cell data");

    int str_to_data = find_data_in_file(filename, "CELL_DATA");

    if (str_to_data == -1)
        write_new_data(filename, "CELL_DATA", dataname, _cells.size(), cell_data);
    else
        write_not_new_data(filename, str_to_data, dataname, _cells.size(), cell_data);
}

NonstatGridSaver::NonstatGridSaver(const GraphGrid& grid, const std::vector<Point2>& nodes_coo, std::string filename)
    : _vtk(GridSaver(grid, nodes_coo)), _file_name(filename), _series_name(_file_name + ".vtk.series")
{
    add_in_series();
}

void NonstatGridSaver::new_time_step(double t)
{
    std::ostringstream oss;
    std::ostringstream fn;
    _times.push_back(t);

    fn << std::setfill('0') << std::setw(4) << _files_count << ".vtk";
    _files_count++;
    _vtk.save_area(fn.str());

    _files.push_back(fn.str());
    add_in_series();
}

void NonstatGridSaver::save_vtk_point_data(const std::vector<double>& data, std::string data_name) const
{
    if (_files.size() > 0)
    {
        _vtk.save_vtk_point_data(data, data_name, _files.back());
    }
}

void NonstatGridSaver::save_vtk_cell_data(const std::vector<double>& data, std::string data_name) const
{
    if (_files.size() > 0)
    {
        _vtk.save_vtk_cell_data(data, data_name, _files.back());
    }
}

std::string NonstatGridSaver::print_files() const
{
    std::string ofs;
    if (_files.size() > 0)
    {
        for (int i = 0; i < _files.size() - 1; i++)
            ofs += "    {\"name\": \"" + _files[i] + "\", \"time\": " + std::to_string(_times[i]) + "}" + ",\n";
        ofs += "    {\"name\": \"" + _files.back() + "\", \"time\": " + std::to_string(_times.back()) + "}\n";
    }
    return ofs;
}

void NonstatGridSaver::add_in_series() const
{
    std::ofstream ofs(_series_name);
    if (!ofs)
        throw std::runtime_error("Failed to open " + _file_name + " for writing");

    ofs << "{" << std::endl;
    ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
    ofs << "  \"files\" : [" << std::endl;
    ofs << print_files();
    ofs << "  ]" << std::endl;
    ofs << "}" << std::endl;
}
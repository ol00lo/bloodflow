#include "tests/test_utils.hpp"

int string_count(std::string file_name)
{
    std::ifstream file(file_name);
    int count = 0;
    std::string line;
    while (std::getline(file, line))
    {
        ++count;
    }
    return count;
}

std::vector<double> result_data(const std::vector<bflow::Point2>& points_coo)
{
    std::vector<double> res_dat(points_coo.size());
    for (int i = 0; i < points_coo.size(); ++i)
    {
        bflow::Point2 p = points_coo[i];
        double r = sqrt(p.x * p.x + p.y * p.y);
        res_dat[i] = cos(r / 5 * 3.1415);
    }
    return res_dat;
}

std::vector<double> distance_from_zero(const bflow::GraphGrid& grid)
{
    std::vector<double> res(grid.n_points());
    double len = 0.0;
    res[0] = len;
    for (int i = 1; i < res.size(); i++)
    {
        len += grid.find_cell_length(i);
        res.push_back(len);
    }
    return res;
}

std::vector<double> find_point_data(const std::vector<bflow::Point2>& points_coo, double t)
{
    std::vector<double> res_dat(points_coo.size());
    for (int i = 0; i < points_coo.size(); ++i)
    {
        bflow::Point2 p = points_coo[i];
        double r = sqrt(p.x * p.x + p.y * p.y);
        res_dat[i] = t*cos(r*t / 5 * 3.1415) ;
    }
    return res_dat;
}
std::vector<double> find_cell_data(const std::vector<double>& cell, double t)
{
    std::vector<double> res_dat(cell.size());
    for (int i = 0; i < cell.size(); ++i)
    {
        res_dat[i] = t*cos(cell[i]*t / 5 * 3.1415);
    }
    return res_dat;
}

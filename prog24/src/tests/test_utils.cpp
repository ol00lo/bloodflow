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
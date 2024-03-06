#ifndef TEST_UTILS
#define TEST_UTILS

#include <fstream>
#include <string>
#include <vector>
#include "bflow/graph_grid.hpp"

int string_count(std::string file_name);
std::vector<double> result_data(const std::vector<bflow::Point2>& points_coo);
std::vector<double> distance_from_zero(const bflow::GraphGrid& grid);
std::vector<double> find_point_data(const std::vector<bflow::Point2>& points_coo, double t);
std::vector<double> find_cell_data(const std::vector<double>& points_coo, double t);

#endif
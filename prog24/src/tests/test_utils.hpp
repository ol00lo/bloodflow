#ifndef TEST_UTILS
#define TEST_UTILS

#include <fstream>
#include <string>
#include <vector>
#include "bflow/graph_grid.hpp"

int string_count(std::string file_name);
std::vector<double> result_data(const std::vector<bflow::Point2>& points_coo);

#endif
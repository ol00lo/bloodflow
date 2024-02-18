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
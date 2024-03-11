#ifndef DEBUG_HPP
#define DEBUG_HPP

#include "bflow/matrix.hpp"
#include <fstream>
#include <iostream>

void print_matrix_full(const bflow::ISparseMatrix& mat, std::ostream& s = std::cout);
void print_matrix_stencil(size_t irow, const bflow::ISparseMatrix& mat, std::ostream& s = std::cout);
void print_matrix_stencil(const bflow::ISparseMatrix& mat, std::ostream& s = std::cout);

#endif
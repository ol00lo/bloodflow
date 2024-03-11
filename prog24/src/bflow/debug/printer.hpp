#ifndef PRINTER_HPP
#define PRINTER_HPP

#include "bflow/matrix.hpp"
#include <iostream>
namespace bflow
{
namespace dbg
{
void print_matrix_full(const ISparseMatrix& mat, std::ostream& s = std::cout);
void print_matrix_stencil(int irow, const ISparseMatrix& mat, std::ostream& s = std::cout);
void print_matrix_stencil(const ISparseMatrix& mat, std::ostream& s = std::cout);
}
}

#endif
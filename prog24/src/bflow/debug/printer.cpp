#include "bflow/debug/printer.hpp"
#include <iomanip>

using namespace bflow;

void dbg::print_matrix_full(const ISparseMatrix& mat, std::ostream& s)
{
    for (int i = 0; i < mat.n_rows(); i++)
    {
        for (int j = 0; j < mat.n_rows(); j++)
        {
            if (mat.is_in_stencil(i, j) == false)
            {
                s << std::setw(4) << "*";
            }
            else
            {
                s << std::setw(4) << mat.value(i, j);
            }
        }
        s << std::endl;
    }
}

void dbg::print_matrix_stencil(int irow, const ISparseMatrix& mat, std::ostream& s)
{
    s << "ROW = " << irow << std::endl;
    for (int j = 0; j < mat.n_rows() - 1; j++)
    {
        if (mat.is_in_stencil(irow, j) == true)
        {
            s << j << ": " << mat.value(irow, j) << std::endl;
        }
    }
}

void dbg::print_matrix_stencil(const ISparseMatrix& mat, std::ostream& s)
{
    s << "NROWS = " << mat.n_rows() << std::endl;
    for (int i = 0; i < mat.n_rows(); i++)
    {
        print_matrix_stencil(i, mat, s);
    }
}
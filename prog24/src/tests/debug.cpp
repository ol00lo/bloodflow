#include "debug.hpp"
#include <iomanip>

void print_matrix_full(const bflow::ISparseMatrix& mat, std::ostream& s)
{
    for (size_t i = 0; i < mat.n_rows(); i++)
    {
        for (size_t j = 0; j < mat.n_rows(); j++)
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

void print_matrix_stencil(size_t irow, const bflow::ISparseMatrix& mat, std::ostream& s)
{
    s << "ROW = " << irow << std::endl;
    for (size_t j = 0; j < mat.n_rows() - 1; j++)
    {
        if (mat.is_in_stencil(irow, j) == true)
        {
            s << j << ": " << mat.value(irow, j) << std::endl;
        }
    }
}

void print_matrix_stencil(const bflow::ISparseMatrix& mat, std::ostream& s)
{
    s << "NROWS = " << mat.n_rows() << std::endl;
    for (size_t i = 0; i < mat.n_rows(); i++)
    {
        print_matrix_stencil(i, mat, s);
    }
}
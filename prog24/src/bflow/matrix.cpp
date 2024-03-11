#include "matrix.hpp"
#include <iostream>

using namespace bflow;
constexpr size_t INVALID_INDEX = size_t(-1);

void ISparseMatrix::validate_ij(size_t i, size_t j) const
{
    if (i >= n_rows() || j >= n_rows())
    {
        throw std::runtime_error("values out of matrix");
    }
}

size_t CsrMatrix::n_rows() const
{
    return _addr.size() - 1;
}

double CsrMatrix::value(size_t i, size_t j) const
{
    size_t a = find_index(i, j);
    if (a != INVALID_INDEX)
        return _vals[a];
    else
        return 0.0;
}

bool CsrMatrix::is_in_stencil(size_t i, size_t j) const
{
    size_t a = find_index(i, j);
    if (a != INVALID_INDEX)
        return true;
    else
        return false;
}

size_t CsrMatrix::find_index(size_t i, size_t j) const
{
    validate_ij(i, j);
    size_t addr0 = _addr[i];
    size_t addr1 = _addr[i + 1];

    for (size_t a = addr0; a < addr1; a++)
    {
        if (_cols[a] == j)
        {
            return a;
        }
    }
    return INVALID_INDEX;
}

void LodMatrix::set_value(size_t i, size_t j, double val)
{
    validate_ij(i, j);
    _data[i][j] = val;
}

void LodMatrix::clear_row(size_t i)
{
    validate_ij(i, 0);
    _data[i].clear();
}

CsrMatrix LodMatrix::to_csr() const
{
    std::vector<size_t> addr;
    std::vector<size_t> cols;
    std::vector<double> vals;
    addr.push_back(0.0);
    int cout = 0;
    for (size_t i = 0; i < _data.size(); ++i)
    {
        for (auto it = _data[i].cbegin(); it != _data[i].cend(); ++it)
        {
            vals.push_back(it->second);
            cols.push_back(it->first);
            cout++;
        }
        addr.push_back(cout);
    }
    CsrMatrix s = CsrMatrix(addr, cols, vals);
    return s;
}

size_t LodMatrix::n_rows() const
{
    return _data.size();
}

double LodMatrix::value(size_t i, size_t j) const
{
    validate_ij(i, j);
    auto fnd = _data[i].find(j);
    if (fnd != _data[i].end())
    {
        return fnd->second;
    }
    else
        return 0.0;
}

bool LodMatrix::is_in_stencil(size_t i, size_t j) const
{
    validate_ij(i, j);
    auto fnd = _data[i].find(j);
    if (fnd != _data[i].end())
        return true;
    else
        return false;
}
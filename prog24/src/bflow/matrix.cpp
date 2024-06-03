#include "matrix.hpp"
#include <iostream>

using namespace bflow;

void ISparseMatrix::validate_ij(int i, int j) const
{
    if (i >= n_rows() || j >= n_rows() || i < 0 || j < 0)
    {
        throw std::runtime_error("values out of matrix");
    }
}

int CsrMatrix::n_rows() const
{
    return _addr.size() - 1;
}

double CsrMatrix::value(int i, int j) const
{
    int a = find_index(i, j);
    if (a != INVALID_INDEX)
        return _vals[a];
    else
        return 0.0;
}

const std::vector<int>& CsrMatrix::addr() const
{
    return _addr;
}

const std::vector<int>& CsrMatrix::cols() const
{
    return _cols;
}

int CsrMatrix::n_nonzeros() const
{
    return _cols.size();
}
std::vector<double>& CsrMatrix::vals()
{
    return _vals;
}
const std::vector<double>& CsrMatrix::vals() const
{
    return _vals;
}

bool CsrMatrix::is_in_stencil(int i, int j) const
{
    int a = find_index(i, j);
    if (a != INVALID_INDEX)
        return true;
    else
        return false;
}

int CsrMatrix::find_index(int i, int j) const
{
    validate_ij(i, j);
    int addr0 = _addr[i];
    int addr1 = _addr[i + 1];

    for (int a = addr0; a < addr1; a++)
    {
        if (_cols[a] == j)
        {
            return a;
        }
    }
    return INVALID_INDEX;
}
void CsrMatrix::set_unit_row(size_t irow)
{
    const std::vector<int>& a = addr();
    const std::vector<int>& c = cols();

    int start = a.at(irow);
    int end = a.at(irow + 1);
    for (int i = start; i < end; ++i)
    {
        _vals[i] = (c[i] == irow) ? 1.0 : 0.0;
    }
}

void CsrMatrix::set_stencil(const std::vector<std::set<int>>& stencil_set)
{
    _addr = std::vector<int>(1, 0);
    _cols.clear();
    for (int irow = 0; irow < stencil_set.size(); ++irow)
    {
        const std::set<int>& cols = stencil_set[irow];
        _addr.push_back(_addr.back() + cols.size());
        for (int col : cols)
        {
            _cols.push_back(col);
            _vals.push_back(0);
        }
    }
}

std::vector<double> CsrMatrix::mult_vec(const std::vector<double>& u) const
{
    const std::vector<int>& a = addr();
    const std::vector<int>& c = cols();
    const std::vector<double>& v = vals();

    std::vector<double> ret(n_rows(), 0);
    for (size_t irow = 0; irow < n_rows(); ++irow)
    {
        size_t start = a[irow];
        size_t end = a[irow + 1];
        for (size_t i = start; i < end; ++i)
        {
            ret[irow] += v[i] * u[c[i]];
        }
    }

    return ret;
}
void LodMatrix::set_value(int i, int j, double val)
{
    validate_ij(i, j);
    _data[i][j] = val;
}

void LodMatrix::clear_row(int i)
{
    validate_ij(i, 0);
    _data[i].clear();
}

CsrMatrix LodMatrix::to_csr() const
{
    std::vector<int> addr;
    std::vector<int> cols;
    std::vector<double> vals;
    addr.push_back(0.0);
    int cout = 0;
    for (int i = 0; i < _data.size(); ++i)
    {
        for (auto it = _data[i].cbegin(); it != _data[i].cend(); ++it)
        {
            vals.push_back(it->second);
            cols.push_back(it->first);
            cout++;
        }
        addr.push_back(cout);
    }
    CsrMatrix s(std::move(addr), std::move(cols), std::move(vals));
    return s;
}

int LodMatrix::n_rows() const
{
    return _data.size();
}

double LodMatrix::value(int i, int j) const
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

bool LodMatrix::is_in_stencil(int i, int j) const
{
    validate_ij(i, j);
    auto fnd = _data[i].find(j);
    if (fnd != _data[i].end())
        return true;
    else
        return false;
}
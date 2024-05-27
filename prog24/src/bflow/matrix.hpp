#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>
#include <set>

namespace bflow
{
constexpr int INVALID_INDEX = int(-1);
class ISparseMatrix
{
public:
    virtual ~ISparseMatrix() = default;

    virtual int n_rows() const = 0;
    virtual double value(int i, int j) const = 0;
    virtual bool is_in_stencil(int i, int j) const = 0;

protected:
    void validate_ij(int i, int j) const;
};

class CsrMatrix : public ISparseMatrix
{
public:
    CsrMatrix() = default;
    CsrMatrix(const std::vector<int>& addr, const std::vector<int>& cols, const std::vector<double>& vals)
        : _addr(addr), _cols(cols), _vals(vals) {}
    CsrMatrix(std::vector<int>&& addr, std::vector<int>&& cols, std::vector<double>&& vals)
        : _addr(std::move(addr)), _cols(std::move(cols)), _vals(std::move(vals)) {}
    int n_rows() const override;
    double value(int i, int j) const override;
    bool is_in_stencil(int i, int j) const override;
    const std::vector<double>& vals() const;
    std::vector<double>& vals();
    int n_nonzeros() const;
    const std::vector<int>& addr() const;
    const std::vector<int>& cols() const;
    int find_index(int i, int j) const;
    void set_unit_row(size_t irow);
    void set_stencil(const std::vector<std::set<int>>& stencil_set);
    std::vector<double> mult_vec(const std::vector<double>& u) const;

private:
    std::vector<int> _addr;
    std::vector<int> _cols;
    std::vector<double> _vals;
};

class LodMatrix : public ISparseMatrix
{
public:
    LodMatrix(int nrows) : _data(nrows)
    {
    }

    void set_value(int i, int j, double val);
    void clear_row(int i);
    CsrMatrix to_csr() const;

    int n_rows() const override;
    double value(int i, int j) const override;
    bool is_in_stencil(int i, int j) const override;

private:
    std::vector<std::map<int, double>> _data;
};
} // namespace bflow

#endif
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>

namespace bflow
{
class ISparseMatrix
{
public:
    virtual ~ISparseMatrix() = default;

    virtual size_t n_rows() const = 0;
    virtual double value(size_t i, size_t j) const = 0;
    virtual bool is_in_stencil(size_t i, size_t j) const = 0;

protected:
    void validate_ij(size_t i, size_t j) const;
};

class CsrMatrix : public ISparseMatrix
{
public:
    CsrMatrix(const std::vector<size_t>& addr, const std::vector<size_t>& cols, const std::vector<double>& vals)
        : _addr(addr), _cols(cols), _vals(vals){}

    size_t n_rows() const override;
    double value(size_t i, size_t j) const override;
    bool is_in_stencil(size_t i, size_t j) const override;

private:
    std::vector<size_t> _addr;
    std::vector<size_t> _cols;
    std::vector<double> _vals;
    size_t find_index(size_t i, size_t j) const;
};

class LodMatrix : public ISparseMatrix
{
public:
    LodMatrix(size_t nrows) : _data(nrows){}

    void set_value(size_t i, size_t j, double val);
    void clear_row(size_t i);
    CsrMatrix to_csr() const;

    size_t n_rows() const override;
    double value(size_t i, size_t j) const override;
    bool is_in_stencil(size_t i, size_t j) const override;

private:
    std::vector<std::map<size_t, double>> _data;
};
} // namespace bflow

#endif
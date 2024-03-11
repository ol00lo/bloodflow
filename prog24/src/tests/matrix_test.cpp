#include "catch.hpp"
#include "tests/debug.hpp"
#include <sstream>

using namespace bflow;

TEST_CASE("test sparse matrix 1", "[sparse_matrix1]")
{
    bool exception_thrown;

    // =========== constructor
    std::vector<size_t> addr{0, 1, 3, 4};
    std::vector<size_t> cols{2, 0, 2, 1};
    std::vector<double> vals{10, 4, 2, 3};
    CsrMatrix cmat(addr, cols, vals);
    // ============ n-rows
    CHECK(cmat.n_rows() == 3);
    // ============ value
    CHECK(cmat.value(0, 0) == 0);
    CHECK(cmat.value(0, 1) == 0);
    CHECK(cmat.value(0, 2) == 10.0);
    CHECK(cmat.value(1, 0) == 4);
    CHECK(cmat.value(1, 1) == 0);
    CHECK(cmat.value(1, 2) == 2);
    CHECK(cmat.value(2, 0) == 0);
    CHECK(cmat.value(2, 1) == 3);
    CHECK(cmat.value(2, 2) == 0);
    try
    {
        exception_thrown = false;
        cmat.value(2, 3);
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    // ============ is_in_stencil
    CHECK(cmat.is_in_stencil(0, 0) == false);
    CHECK(cmat.is_in_stencil(0, 2) == true);
    CHECK(cmat.is_in_stencil(2, 2) == false);
    try
    {
        exception_thrown = false;
        cmat.value(3, 0);
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);

    // LODMatrix
    LodMatrix lmat(3);
    // nrows
    CHECK(lmat.n_rows() == 3);
    CHECK(lmat.value(0, 2) == 0);
    // set_value
    lmat.set_value(0, 2, 10.0);
    lmat.set_value(1, 0, 4.0);
    lmat.set_value(1, 2, 2.0);
    lmat.set_value(2, 1, 3.0);
    try
    {
        exception_thrown = false;
        lmat.set_value(10, 10, 1.0);
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    // value
    CHECK(lmat.value(0, 2) == 10.0);
    CHECK(lmat.value(2, 1) == 3.0);
    CHECK(lmat.value(2, 2) == 0.0);
    try
    {
        exception_thrown = false;
        lmat.value(1, 3);
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    // convert
    CsrMatrix cmat2 = lmat.to_csr();
    // check conversion through print
    std::ostringstream ss_lmat, ss_cmat, ss_cmat2;
    print_matrix_full(lmat, ss_lmat);
    print_matrix_full(cmat, ss_cmat);
    print_matrix_full(cmat2, ss_cmat2);
    CHECK(ss_lmat.str() == ss_cmat2.str());
    CHECK(ss_cmat.str() == ss_cmat2.str());
}

TEST_CASE("test sparse matrix 2", "[sparse_matrix2]")
{
    bool exception_thrown;

    LodMatrix lmat(3);
    lmat.set_value(0, 0, 1);
    lmat.set_value(1, 0, 2);
    lmat.set_value(2, 0, 3);
    try
    {
        lmat.set_value(3, 0, 4);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    try
    {
        lmat.set_value(0, 3, 4);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    lmat.clear_row(2);
    try
    {
        lmat.clear_row(3);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    CHECK(lmat.value(0, 0) == 1);
    CHECK(lmat.value(1, 0) == 2);
    CHECK(lmat.is_in_stencil(2, 0) == false);
    CHECK(lmat.value(2, 0) == 0);

    lmat.set_value(0, 0, 0);
    CHECK(lmat.value(0, 0) == 0);

    CsrMatrix cmat = lmat.to_csr();
    // CsrMatrix cmat({0, 1, 1, 2}, {0, 0}, {1, 2});

    CHECK(cmat.value(0, 0) == 0);
    CHECK(cmat.value(1, 0) == 2);
    CHECK(cmat.is_in_stencil(2, 0) == false);
    CHECK(cmat.value(2, 0) == 0);
    try
    {
        cmat.value(3, 0);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    try
    {
        cmat.is_in_stencil(0, 3);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);

    CHECK(cmat.is_in_stencil(1, 2) == false);
    CHECK(cmat.value(1, 2) == 0);

    std::ostringstream oss_cmat, oss_lmat;
    print_matrix_full(lmat, oss_lmat);
    print_matrix_full(cmat, oss_cmat);
    CHECK(oss_lmat.str() == oss_cmat.str());
}

TEST_CASE("test sparse matrix 3", "[sparse_matrix3]")
{
    bool exception_thrown;

    LodMatrix lmat(0);
    CHECK(lmat.n_rows() == 0);
    try
    {
        lmat.value(0, 0);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    try
    {
        lmat.set_value(0, 0, 0);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    try
    {
        lmat.clear_row(0);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);

    CsrMatrix cmat({0}, {}, {});
    CHECK(cmat.n_rows() == 0);
    try
    {
        cmat.value(0, 0);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
    try
    {
        cmat.is_in_stencil(0, 0);
        exception_thrown = false;
    }
    catch (std::runtime_error&)
    {
        exception_thrown = true;
    }
    CHECK(exception_thrown);
}
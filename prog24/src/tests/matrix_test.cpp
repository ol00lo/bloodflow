#include "catch.hpp"
#include "bflow/debug/printer.hpp"
#include <sstream>

using namespace bflow;

TEST_CASE("test sparse matrix 1", "[sparse_matrix1]")
{
    // =========== constructor
    std::vector<int> addr{0, 1, 3, 4};
    std::vector<int> cols{2, 0, 2, 1};
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
    CHECK_THROWS(cmat.value(2, 3));
    // ============ is_in_stencil
    CHECK(cmat.is_in_stencil(0, 0) == false);
    CHECK(cmat.is_in_stencil(0, 2) == true);
    CHECK(cmat.is_in_stencil(2, 2) == false);
    CHECK_THROWS(cmat.value(3, 0));

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
    CHECK_THROWS(lmat.set_value(10, 10, 1.0));
    // value
    CHECK(lmat.value(0, 2) == 10.0);
    CHECK(lmat.value(2, 1) == 3.0);
    CHECK(lmat.value(2, 2) == 0.0);
    CHECK_THROWS(lmat.value(1, 3));
    // convert
    CsrMatrix cmat2 = lmat.to_csr();
    // check conversion through print
    std::ostringstream ss_lmat, ss_cmat, ss_cmat2;
    dbg::print_matrix_full(lmat, ss_lmat);
    dbg::print_matrix_full(cmat, ss_cmat);
    dbg::print_matrix_full(cmat2, ss_cmat2);
    CHECK(ss_lmat.str() == ss_cmat2.str());
    CHECK(ss_cmat.str() == ss_cmat2.str());
}

TEST_CASE("test sparse matrix 2", "[sparse_matrix2]")
{
    LodMatrix lmat(3);
    lmat.set_value(0, 0, 1);
    lmat.set_value(1, 0, 2);
    lmat.set_value(2, 0, 3);
    CHECK_THROWS(lmat.set_value(3, 0, 4));
    CHECK_THROWS(lmat.set_value(0, 3, 4));

    lmat.clear_row(2);
    CHECK_THROWS(lmat.clear_row(3));

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
    CHECK_THROWS(cmat.value(3, 0));
    CHECK_THROWS(cmat.is_in_stencil(0, 3));

    CHECK(cmat.is_in_stencil(1, 2) == false);
    CHECK(cmat.value(1, 2) == 0);

    std::ostringstream oss_cmat, oss_lmat;
    dbg::print_matrix_full(lmat, oss_lmat);
    dbg::print_matrix_full(cmat, oss_cmat);
    CHECK(oss_lmat.str() == oss_cmat.str());
}

TEST_CASE("test sparse matrix 3", "[sparse_matrix3]")
{
    LodMatrix lmat(0);
    CHECK(lmat.n_rows() == 0);
    CHECK_THROWS(lmat.value(0, 0));
    CHECK_THROWS(lmat.set_value(0, 0, 0));
    CHECK_THROWS(lmat.clear_row(0));

    CsrMatrix cmat({0}, {}, {});
    CHECK(cmat.n_rows() == 0);
    CHECK_THROWS(cmat.value(0,0));
    CHECK_THROWS(cmat.is_in_stencil(0, 0));
}
#include "bflow/graph_grid.hpp"
#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/vessel_graph.hpp"
#include "bflow/vtk.hpp"
#include "tests/nonlinsolver.hpp"
#include "catch.hpp"
#include <fstream>

#include <iomanip>
using namespace bflow;
namespace
{

class FemGrid
{
public:
    FemGrid(double length, size_t n, size_t power) : _power(power)
    {
        _points.push_back(0);
        _nodes.push_back(0);
        for (size_t i = 0; i < n; ++i)
        {
            _points.push_back(double(i + 1) / n * length);
            _nodes.push_back(_points.back());
            _nodes.push_back(_points.back());
        }
        _nodes.pop_back();
        _h = _points.back() / n_elements();
        if (power > 1)
        {
            _THROW_NOT_IMP_;
        }
    }
    FemGrid(const GraphGrid& grid, std::vector<Point2> nodes_coo) : _power(grid.n_midnodes)
    {
        if (_power > 6)
        {
            throw std::runtime_error("Only power=1,2,3,4,5,6 is allowed");
        }
        if (grid.n_edges() != 1)
        {
            throw std::runtime_error("Only single vessel grids are allowed");
        }
        _points.resize(grid.n_points(), 0.0);
        double x = 0;
        for (size_t i = 0; i < grid.n_elem(); ++i)
        {
            x += grid.find_cell_length(i);
            _points[i + 1] = x;
        }
        _nodes.resize(grid.n_nodes() - (_power - 1) * (_points.size() - 1), 0);
        _h = grid.find_cell_length(0);

        for (size_t inode = 0; inode < _nodes.size(); ++inode)
        {
            size_t ipoint = (inode + 1) / 2;
            _nodes[inode] = _points[ipoint];
        }
        double l = _h / _power;
        for (size_t i = 0; i < _points.size() - 1; i++)
        {
            double p = _points[i] + l;
            for (size_t j = 1; j < _power; j++)
            {
                _nodes.push_back(p);
                p += l;
            }
        }
        fill_f_vec();
    }

    double h() const
    {
        return _h;
    }

    size_t n_nodes() const
    {
        return _nodes.size();
    }

    size_t n_elements() const
    {
        return _points.size() - 1;
    }

    size_t n_points() const
    {
        return _points.size();
    }

    size_t n_local_bases() const
    {
        return _power + 1;
    }

    double node(size_t i) const
    {
        return _nodes[i];
    }

    CsrMatrix mass_matrix() const
    {
        CsrMatrix ret(stencil());
        std::vector<double> local = local_mass_matrix();

        for (size_t ielem = 0; ielem < n_elements(); ++ielem)
        {
            std::vector<int> lg = tab_elem_nodes(ielem);

            for (size_t irow = 0; irow < n_local_bases(); ++irow)
                for (size_t icol = 0; icol < n_local_bases(); ++icol)
                {
                    double v = h() / 2 * local[irow * n_local_bases() + icol];
                    size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                    ret.vals()[iaddr] += v;
                }
        }

        return ret;
    }

    CsrMatrix transport_matrix() const
    {
        CsrMatrix ret(stencil());
        std::vector<double> local = local_transport_matrix();

        for (size_t ielem = 0; ielem < n_elements(); ++ielem)
        {
            std::vector<int> lg = tab_elem_nodes(ielem);

            // block diagonal
            for (size_t irow = 0; irow < n_local_bases(); ++irow)
                for (size_t icol = 0; icol < n_local_bases(); ++icol)
                {
                    double v = local[irow * n_local_bases() + icol];
                    size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                    ret.vals()[iaddr] -= v;
                }

            // upwind coupling
            if (ielem > 0)
            {
                // left
                size_t iaddr = ret.find_index(lg[0], lg[0] - 1);
                ret.vals()[iaddr] -= 1;
            }
            {
                // right
                size_t iaddr = ret.find_index(lg[1], lg[1]);
                ret.vals()[iaddr] += 1;
            }
        }

        return ret;
    }

    std::vector<double> load_vector() const
    {
        return _f_vec;
    }

    std::vector<int> tab_elem_nodes(int ielem) const
    {
        std::vector<int> ret;
        if (_power == 1)
        {
            ret = {2 * ielem, 2 * ielem + 1};
        }
        else if (_power == 2)
        {
            ret = {2 * ielem, 2 * ielem + 1, 2 * int(n_elements()) + ielem};
        }
        else if (_power == 3)
        {
            ret = {2 * ielem, 2 * ielem + 1, 2 * int(n_elements()) + 2 * ielem, 2 * int(n_elements()) + 2 * ielem + 1};
        }
        else if (_power == 4)
        {
            ret = {2 * ielem, 2 * ielem + 1, 2 * int(n_elements()) + 3 * ielem, 2 * int(n_elements()) + 3 * ielem + 1,
                   2 * int(n_elements()) + 3 * ielem + 2};
        }
        else if (_power == 5)
        {
            ret = {2 * ielem,
                   2 * ielem + 1,
                   2 * int(n_elements()) + 4 * ielem,
                   2 * int(n_elements()) + 4 * ielem + 1,
                   2 * int(n_elements()) + 4 * ielem + 2,
                   2 * int(n_elements()) + 4 * ielem + 3};
        }
        else if (_power == 6)
        {
            ret = {2 * ielem,
                   2 * ielem + 1,
                   2 * int(n_elements()) + 5 * ielem,
                   2 * int(n_elements()) + 5 * ielem + 1,
                   2 * int(n_elements()) + 5 * ielem + 2,
                   2 * int(n_elements()) + 5 * ielem + 3,
                   2 * int(n_elements()) + 5 * ielem + 4};
        }
        else
        {
            _THROW_NOT_IMP_;
        }
        return ret;
    }

private:
    const int _power;
    mutable CsrMatrix _stencil;
    std::vector<double> _points;
    std::vector<double> _nodes;
    std::vector<double> _f_vec;
    double _h;

    CsrMatrix stencil() const
    {
        if (_stencil.n_rows() < 1)
        {
            std::vector<std::set<int>> lod(n_nodes());
            // block diagonal
            for (size_t ielem = 0; ielem < n_elements(); ++ielem)
            {
                std::vector<int> lg = tab_elem_nodes(ielem);

                for (size_t irow = 0; irow < n_local_bases(); ++irow)
                    for (size_t icol = 0; icol < n_local_bases(); ++icol)
                    {
                        lod[lg[irow]].insert(lg[icol]);
                    }
            }
            // coupling
            for (int ipoint = 0; ipoint < n_points() - 1; ++ipoint)
            {
                std::vector<int> nodes = tab_point_nodes(ipoint);
                if (nodes.size() == 2)
                {
                    lod[nodes[0]].insert(nodes[1]);
                    lod[nodes[1]].insert(nodes[0]);
                }
            }
            _stencil.set_stencil(lod);
        }
        return _stencil;
    }

    std::vector<int> tab_point_nodes(int ipoint) const
    {
        if (ipoint == 0)
            return {0};
        else if (ipoint == n_points() - 1)
            return {2 * ipoint};
        else
            return {2 * ipoint - 1, 2 * ipoint};
    }

    std::vector<double> local_mass_matrix() const
    {
        if (_power == 1)
        {
            return {2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0};
        }
        else if (_power == 2)
        {
            return {4.0 / 15.0, -1.0 / 15.0, 2.0 / 15.0, -1.0 / 15.0, 4.0 / 15.0,
                    2.0 / 15.0, 2.0 / 15.0,  2.0 / 15.0, 16.0 / 15.0};
        }
        else if (_power == 3)
        {
            return {16.0 / 105.0, 19.0 / 840.0, 33.0 / 280.0,  -3.0 / 70,   19.0 / 840.0, 16.0 / 105.0,
                    -3.0 / 70.0,  33.0 / 280.0, 33.0 / 280.0,  -3.0 / 70.0, 27.0 / 35.0,  -27.0 / 280.0,
                    -3.0 / 70.0,  33.0 / 280.0, -27.0 / 280.0, 27.0 / 35.0};
        }
        else if (_power == 4)
        {
            return {292.0 / 2835, -29.0 / 2835, 296.0 / 2835, -58.0 / 945,  8.0 / 405,   -29.0 / 2835, 292.0 / 2835,
                    8.0 / 405,    -58.0 / 945,  296.0 / 2835, 296.0 / 2835, 8.0 / 405,   256.0 / 405,  -128.0 / 945,
                    256.0 / 2835, -58.0 / 945,  -58.0 / 945,  -128.0 / 945, 208.0 / 315, -128.0 / 945, 8.0 / 405,
                    296.0 / 2835, 256.0 / 2835, -128.0 / 945, 256.0 / 405};
        }
        else if (_power == 5)
        {
            return {1907.0 / 24948,    493.0 / 88704,     24775.0 / 266112, -9925.0 / 133056,  -1525.0 / 133056,
                    17125.0 / 399168,  493.0 / 88704,     1907.0 / 24948,   -1525.0 / 133056,  17125.0 / 399168,
                    24775.0 / 266112,  -9925.0 / 133056,  24775.0 / 266112, -1525.0 / 133056,  111625.0 / 199584,
                    -24625.0 / 133056, -62875.0 / 798336, 2125.0 / 14784,   -9925.0 / 133056,  17125.0 / 399168,
                    -24625.0 / 133056, 62375.0 / 99792,   2125.0 / 14784,   -13625.0 / 66528,  -1525.0 / 133056,
                    24775.0 / 266112,  -62875.0 / 798336, 2125.0 / 14784,   111625.0 / 199584, -24625.0 / 133056,
                    17125.0 / 399168,  -9925.0 / 133056,  2125.0 / 14784,   -13625.0 / 66528,  -24625.0 / 133056,
                    62375.0 / 99792};
        }
        else if (_power == 6)
        {
            return {90269.0 / 1501500, -10237.0 / 3003000, 42087.0 / 500500,   -16971.0 / 200200, 10237.0 / 150150,
                    -687.0 / 20020,    3867.0 / 500500,    -10237.0 / 3003000, 90269.0 / 1501500, 3867.0 / 500500,
                    -687.0 / 20020,    10237.0 / 150150,   -16971.0 / 200200,  42087.0 / 500500,  42087.0 / 500500,
                    3867.0 / 500500,   64692.0 / 125125,   -4887.0 / 20020,    5688.0 / 25025,    -14607.0 / 100100,
                    8532.0 / 125125,   -16971.0 / 200200,  -687.0 / 20020,     -4887.0 / 20020,   2619.0 / 4004,
                    -3231.0 / 10010,   9693.0 / 40040,     -14607.0 / 100100,  10237.0 / 150150,  10237.0 / 150150,
                    5688.0 / 25025,    -3231.0 / 10010,    10544.0 / 15015,    -3231.0 / 10010,   5688.0 / 25025,
                    -687.0 / 20020,    -16971.0 / 200200,  -14607.0 / 100100,  9693.0 / 40040,    -3231.0 / 10010,
                    2619.0 / 4004,     -4887.0 / 20020,    3867.0 / 500500,    42087.0 / 500500,  8532.0 / 125125,
                    -14607.0 / 100100, 5688.0 / 25025,     -4887.0 / 20020,    64692.0 / 125125};
        }
        else
        {
            _THROW_NOT_IMP_;
        }
    }

    std::vector<double> local_transport_matrix() const
    {
        if (_power == 1)
        {
            return {-0.5, -0.5, 0.5, 0.5};
        }
        else if (_power == 2)
        {
            return {-1.0 / 2.0, 1.0 / 6.0, -2.0 / 3.0, -1.0 / 6.0, 1.0 / 2.0, 2.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0, 0.0};
        }
        else if (_power == 3)
        {
            return {-1.0 / 2.0,  -7.0 / 80.0,  -57.0 / 80.0, 3.0 / 10.0, 7.0 / 80.0, 1.0 / 2.0,
                    -3.0 / 10.0, 57.0 / 80.0,  57.0 / 80.0,  3.0 / 10.0, 0.0,        -81.0 / 80.0,
                    -3.0 / 10.0, -57.0 / 80.0, 81.0 / 80.0,  0.0};
        }
        else if (_power == 4)
        {
            return {-1.0 / 2,     107.0 / 1890, -736.0 / 945, 134.0 / 315, -64.0 / 315, -107.0 / 1890, 1.0 / 2,
                    64.0 / 315,   -134.0 / 315, 736.0 / 945,  736.0 / 945, -64.0 / 315, 0.0,           -352.0 / 315,
                    512.0 / 945,  -134.0 / 315, 134.0 / 315,  352.0 / 315, 0.0,         -352.0 / 315,  64.0 / 315,
                    -736.0 / 945, -512.0 / 945, 352.0 / 315,  0.0};
        }
        else if (_power == 5)
        {
            return {-1.0 / 2,
                    -5951.0 / 145152,
                    -123425.0 / 145152,
                    6325.0 / 10368,
                    1675.0 / 10368,
                    -575.0 / 1512,
                    5951.0 / 145152,
                    1.0 / 2,
                    -1675.0 / 10368,
                    575.0 / 1512,
                    123425.0 / 145152,
                    -6325.0 / 10368,
                    123425.0 / 145152,
                    1675.0 / 10368,
                    0.0,
                    -6875.0 / 5184,
                    -19375.0 / 48384,
                    51875.0 / 72576,
                    -6325.0 / 10368,
                    -575.0 / 1512,
                    6875.0 / 5184,
                    0.0,
                    51875.0 / 72576,
                    -38125.0 / 36288,
                    -1675.0 / 10368,
                    -123425.0 / 145152,
                    19375.0 / 48384,
                    -51875.0 / 72576,
                    0.0,
                    6875.0 / 5184,
                    575.0 / 1512,
                    6325.0 / 10368,
                    -51875.0 / 72576,
                    38125.0 / 36288,
                    -6875.0 / 5184,
                    0.0};
        }
        else if (_power == 6)
        {
            return {-1.0 / 2,
                    587.0 / 18480,
                    -1776.0 / 1925,
                    5151.0 / 6160,
                    -3967.0 / 5775,
                    732.0 / 1925,
                    -267.0 / 1925,
                    -587.0 / 18480,
                    1.0 / 2,
                    267.0 / 1925,
                    -732.0 / 1925,
                    3967.0 / 5775,
                    -5151.0 / 6160,
                    1776.0 / 1925,
                    1776.0 / 1925,
                    -267.0 / 1925,
                    0.0,
                    -3078.0 / 1925,
                    2136.0 / 1925,
                    -243.0 / 385,
                    648.0 / 1925,
                    -5151.0 / 6160,
                    732.0 / 1925,
                    3078.0 / 1925,
                    0.0,
                    -87.0 / 77,
                    3807.0 / 6160,
                    -243.0 / 385,
                    3967.0 / 5775,
                    -3967.0 / 5775,
                    -2136.0 / 1925,
                    87.0 / 77,
                    0.0,
                    -87.0 / 77,
                    2136.0 / 1925,
                    -732.0 / 1925,
                    5151.0 / 6160,
                    243.0 / 385,
                    -3807.0 / 6160,
                    87.0 / 77,
                    0.0,
                    -3078.0 / 1925,
                    267.0 / 1925,
                    -1776.0 / 1925,
                    -648.0 / 1925,
                    243.0 / 385,
                    -2136.0 / 1925,
                    3078.0 / 1925,
                    0.0};
        }
        else
        {
            _THROW_NOT_IMP_;
        }
    }

    void fill_f_vec()
    {
        if (_power == 1)
        {
            for (size_t i = 0; i < _nodes.size(); i++)
                _f_vec.push_back(_h / 2);
        }
        if (_power == 2)
        {
            for (size_t i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(_h / 3 / 2);
                _f_vec.push_back(_h / 3 / 2);
            }
            for (int i = 0; i < _points.size(); i++)
                _f_vec.push_back(4.0 * _h / 3 / 2);
        }
        if (_power == 3)
        {
            for (size_t i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(_h / 4 / 2);
                _f_vec.push_back(_h / 4 / 2);
            }
            for (int i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(3.0 * _h / 4 / 2);
                _f_vec.push_back(3.0 * _h / 4 / 2);
            }
        }
        if (_power == 4)
        {
            for (size_t i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(7.0 * _h / 45 / 2);
                _f_vec.push_back(7.0 * _h / 45 / 2);
            }
            for (int i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(32.0 * _h / 45 / 2);
                _f_vec.push_back(4.0 * _h / 15 / 2);
                _f_vec.push_back(32.0 * _h / 45 / 2);
            }
        }
        if (_power == 5)
        {
            for (size_t i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(19.0 * _h / 144 / 2);
                _f_vec.push_back(19.0 * _h / 144 / 2);
            }
            for (int i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(25.0 * _h / 48 / 2);
                _f_vec.push_back(25.0 * _h / 72 / 2);
                _f_vec.push_back(25.0 * _h / 48 / 2);
                _f_vec.push_back(25.0 * _h / 72 / 2);
            }
        }
        if (_power == 6)
        {
            for (size_t i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(41.0 * _h / 420 / 2);
                _f_vec.push_back(41.0 * _h / 420 / 2);
            }
            for (int i = 0; i < _points.size(); i++)
            {
                _f_vec.push_back(18.0 * _h / 35 / 2);
                _f_vec.push_back(9.0 * _h / 140 / 2);
                _f_vec.push_back(68.0 * _h / 105 / 2);
                _f_vec.push_back(9.0 * _h / 140 / 2);
                _f_vec.push_back(18.0 * _h / 35 / 2);
            }
        }
    }
};

struct ElementBoundaryFluxes
{
    double a_x0 = 0;
    double a_x1 = 0;
    double u_x0 = 0;
    double u_x1 = 0;
};

struct ProblemData
{
    ProblemData()
    {
        recompute();
    }

    void recompute()
    {
        beta = 4.0 / 3.0 * sqrt(pi) * h * E / area0;
        amult = 4 * sqrt(beta / 2 / rho);
        root4_a0 = sqrt(sqrt(area0));
        visc_coef = -2 * (profile_order + 2) * mu * pi / profile_order / rho;
    }

    static constexpr double pi = 3.1415926;

    // geometry parameters
    double L = 1;
    double area0 = pi * 1e-4;

    // fluid parameters
    double rho = 1050;
    double mu = 0;
    double profile_order = 9;

    // vessel tissue parameters
    double h = 1.5e-3;
    double E = 4e5;

    // inflow conditions
    double q_inflow(double t) const
    {
        return 1e-6 * exp(-1e4 * (t - 0.05) * (t - 0.05));
    };

    // p(a)
    double pressure(double area) const
    {
        return 0 + beta * (std::sqrt(area) - std::sqrt(area0));
    };

    // fluxes
    double flux_a(double area, double vel) const
    {
        return area * vel;
    };
    double flux_u(double area, double vel) const
    {
        return vel * vel / 2 + pressure(area) / rho;
    }

    // characteristic variables
    double w1(double area, double vel) const
    {
        return vel + amult * (sqrt(sqrt(area)) - root4_a0);
    }
    double w2(double area, double vel) const
    {
        return vel - amult * (sqrt(sqrt(area)) - root4_a0);
    }

    // computed characteristics
    double beta;
    double amult;
    double root4_a0;
    double visc_coef;
};

class IUpwindFluxCalculator
{
public:
    virtual void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                         std::vector<ElementBoundaryFluxes>& fluxes) = 0;
};

class InflowQFluxCalculator : public IUpwindFluxCalculator
{
public:
    InflowQFluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> get_time, size_t cell)
        : _data(data), _get_time(get_time), _cell_right(cell), _node_right(grid.tab_elem_nodes(cell)[0]),
          _sys(data.amult, data.root4_a0)
    {
        _a = data.area0;
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double t = _get_time();
        double q = _data.q_inflow(t);
        double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);
        // double w2_upw = 0;

        double area_upw = _a;
        double velo_upw = q / area_upw;
        _sys.set_qw(q, w2_upw);
        solve_nonlinear_system(_sys, area_upw, velo_upw, 1e-12);

        double fa = _data.flux_a(area_upw, velo_upw);
        // double fa = q;
        double fu = _data.flux_u(area_upw, velo_upw);

        fluxes[_cell_right].a_x0 = fa;
        fluxes[_cell_right].u_x0 = fu;

        _a = area_upw;
    }

private:
    struct NonlinearSystem : public INonlinearSystem2
    {
    public:
        NonlinearSystem(double mult, double root4_area0) : _mult(mult), _root4_area0(root4_area0)
        {
            set_qw(0, 0);
        }

        void set_qw(double q, double w2)
        {
            _q = q;
            _w2 = w2;
        };

        std::array<double, 2> f(double area, double velo) const override
        {
            return {area * velo - _q, velo - _mult * (sqrt(sqrt(area)) - _root4_area0) - _w2};
        };
        std::array<double, 4> jac(double area, double velo) const override
        {
            return {velo, area, -0.25 * _mult * pow(area, -0.75), 1};
        };

    private:
        const double _mult, _root4_area0;
        double _q, _w2;
    };

    const ProblemData& _data;
    std::function<double()> _get_time;
    const size_t _cell_right;
    const size_t _node_right;

    NonlinearSystem _sys;
    double _a;
};

class OutflowFluxCalculator : public IUpwindFluxCalculator
{
public:
    OutflowFluxCalculator(const FemGrid& grid, const ProblemData& data, size_t cell)
        : _data(data), _cell_left(cell), _node_left(grid.tab_elem_nodes(cell)[1])
    {
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double w1_upw = _data.w1(area[_node_left], velocity[_node_left]);
        double w2_upw = 0;

        double a1 = (w1_upw - w2_upw) / 2 / _data.amult + _data.root4_a0;
        double area_upw = a1 * a1 * a1 * a1;
        double velo_upw = (w1_upw + w2_upw) / 2.0;

        double fa = _data.flux_a(area_upw, velo_upw);
        double fu = _data.flux_u(area_upw, velo_upw);

        fluxes[_cell_left].a_x1 = fa;
        fluxes[_cell_left].u_x1 = fu;
    }

private:
    const ProblemData& _data;
    const size_t _cell_left;
    const size_t _node_left;
};

class InternalFluxCalculator : public IUpwindFluxCalculator
{
public:
    InternalFluxCalculator(const FemGrid& grid, const ProblemData& data, size_t cell_left, size_t cell_right)
        : _data(data), _node_left(grid.tab_elem_nodes(cell_left)[1]), _node_right(grid.tab_elem_nodes(cell_right)[0]),
          _cell_left(cell_left), _cell_right(cell_right), _amult(data.rho * data.rho / 4.0 / data.beta / data.beta)
    {
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double w1_upw = _data.w1(area[_node_left], velocity[_node_left]);
        double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);
        double a1 = (w1_upw - w2_upw) / 2 / _data.amult + _data.root4_a0;
        double area_upw = a1 * a1 * a1 * a1;
        double velo_upw = (w1_upw + w2_upw) / 2.0;

        double fa = _data.flux_a(area_upw, velo_upw);
        double fu = _data.flux_u(area_upw, velo_upw);

        fluxes[_cell_left].a_x1 = fa;
        fluxes[_cell_right].a_x0 = fa;
        fluxes[_cell_left].u_x1 = fu;
        fluxes[_cell_right].u_x0 = fu;
    }

private:
    const ProblemData& _data;
    const size_t _node_left;
    const size_t _node_right;
    const size_t _cell_left;
    const size_t _cell_right;
    const double _amult;
};

class MergingFluxCalculator : public IUpwindFluxCalculator
{
public:
    MergingFluxCalculator(const FemGrid& grid, const ProblemData& datal, const ProblemData& datar, size_t cell_left,
                          size_t cell_right)
        : _data1(datal), _data2(datar), _node_left(grid.tab_elem_nodes(cell_left)[1]),
          _node_right(grid.tab_elem_nodes(cell_right)[0]), _cell_left(cell_left), _cell_right(cell_right),
          _sys(4 * sqrt(datal.beta / 2 / datal.rho), 4 * sqrt(datar.beta / 2 / datar.rho), datal.root4_a0,
               datar.root4_a0, datal.beta, datar.beta, sqrt(datal.area0), sqrt(datar.area0))
    {
        _area1_upw = _data1.area0;
        _area2_upw = _data2.area0;
        _velo1_upw = 0;
        _velo2_upw = 0;
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double w1_upw = _data1.w1(area[_node_left], velocity[_node_left]);
        double w2_upw = _data2.w2(area[_node_right], velocity[_node_right]);

        _sys.set_ww(w1_upw, w2_upw);
        solve_nonlinear_system(_sys, _area1_upw, _velo1_upw, _area2_upw, _velo2_upw, 1e-12);

        double fa1 = _data1.flux_a(_area1_upw, _velo1_upw);
        double fu1 = _data1.flux_u(_area1_upw, _velo1_upw);
        double fa2 = _data2.flux_a(_area2_upw, _velo2_upw);
        double fu2 = _data2.flux_u(_area2_upw, _velo2_upw);

        fluxes[_cell_left].a_x1 = fa1;
        fluxes[_cell_right].a_x0 = fa2;
        fluxes[_cell_left].u_x1 = fu1;
        fluxes[_cell_right].u_x0 = fu2;
    }

private:
    struct NonlinearSystem : public INonlinearSystem4
    {
    public:
        NonlinearSystem(double mult1, double mult2, double root4_area0_1, double root4_area0_2, double beta1,
                        double beta2, double root2_area0_1, double root2_area0_2)
            : _mult1(mult1), _root4_area0_1(root4_area0_1), _bet1(beta1), _root2_area0_1(root2_area0_1), _mult2(mult2),
              _root4_area0_2(root4_area0_2), _bet2(beta2), _root2_area0_2(root2_area0_2)
        {
            set_ww(0, 0);
        }

        void set_ww(double w1, double w2)
        {
            _w1 = w1;
            _w2 = w2;
        };

        std::array<double, 4> f(double area1, double velo1, double area2, double velo2) const override
        {
            return {velo1 + _mult1 * (std::pow(area1, 0.25) - _root4_area0_1) - _w1,
                    velo2 - _mult2 * (std::pow(area2, 0.25) - _root4_area0_2) - _w2, area1 * velo1 - area2 * velo2,
                    velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo2 * velo2 / 2 -
                        _bet2 * (sqrt(area2) - _root2_area0_2)};
        };
        std::array<double, 16> jac(double area1, double velo1, double area2, double velo2) const override
        {
            return {std::pow(area1, -0.75) * _mult1 / 4,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -std::pow(area2, -0.75) * _mult2 / 4,
                    1.0,
                    velo1,
                    area1,
                    -velo2,
                    -area2,
                    std::pow(area1, -0.5) * _bet1 / 2,
                    velo1,
                    -std::pow(area2, -0.5) * _bet2 / 2,
                    -velo2};
            // return{
            //	std::pow(area1,-0.75) * _mult1 / 4, 0.0, velo1, std::pow(area1,-0.5) * _bet1 / 2,
            //	1.0, 0.0, area1, velo1,
            //	-std::pow(area2,-0.75) * _mult2 / 4, 0.0, -velo2, -std::pow(area2,-0.5) * _bet2 / 2,
            //	1.0, 0.0, -area2, velo2
            //};
        };

    private:
        double _mult1, _root4_area0_1, _bet1, _root2_area0_1;
        double _mult2, _root4_area0_2, _bet2, _root2_area0_2;
        double _w1, _w2;
    };
    const ProblemData& _data1;
    const ProblemData& _data2;
    const size_t _node_left;
    const size_t _node_right;
    const size_t _cell_left;
    const size_t _cell_right;

    NonlinearSystem _sys;

    double _area1_upw;
    double _area2_upw;
    double _velo1_upw;
    double _velo2_upw;
};

class Bifurcation3FluxCalculator : public IUpwindFluxCalculator
{
    Bifurcation3FluxCalculator(const FemGrid& grid, const ProblemData& data1, const ProblemData& data2,
                               const ProblemData& data3, size_t cell_left, size_t cell_right1, size_t cell_right2)
        : _data1(data1), _data2(data2), _data3(data3), _node_left(grid.tab_elem_nodes(cell_left)[1]),
          _node_right1(grid.tab_elem_nodes(cell_right1)[0]), _node_right2(grid.tab_elem_nodes(cell_right2)[0]),
          _cell_left(cell_left), _cell_right1(cell_right1), _cell_right2(cell_right2),
          _sys(4 * sqrt(data1.beta / 2 / data1.rho), 4 * sqrt(data1.beta / 2 / data1.rho),
               4 * sqrt(data3.beta / 2 / data3.rho), data1.root4_a0, data2.root4_a0, data3.root4_a0, data1.beta,
               data2.beta, data3.beta, sqrt(data1.area0), sqrt(data2.area0), sqrt(data3.area0))

    {
        _area1_upw = _data1.area0;
        _area2_upw = _data2.area0;
        _area3_upw = _data3.area0;
        _velo1_upw = 0;
        _velo2_upw = 0;
        _velo3_upw = 0;
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double w1_upw = _data1.w1(area[_node_left], velocity[_node_left]);
        double w2_upw = _data2.w2(area[_node_right1], velocity[_node_right1]);
        double w3_upw = _data2.w2(area[_node_right2], velocity[_node_right2]);

        _sys.set_ww(w1_upw, w2_upw, w3_upw);
        solve_nonlinear_system(_sys, _area1_upw, _velo1_upw, _area2_upw, _velo2_upw, _area3_upw, _velo3_upw, 1e-12);

        double fa1 = _data1.flux_a(_area1_upw, _velo1_upw);
        double fu1 = _data1.flux_u(_area1_upw, _velo1_upw);
        double fa2 = _data2.flux_a(_area2_upw, _velo2_upw);
        double fu2 = _data2.flux_u(_area2_upw, _velo2_upw);
        double fa3 = _data2.flux_a(_area3_upw, _velo3_upw);
        double fu3 = _data2.flux_u(_area3_upw, _velo3_upw);

        fluxes[_cell_left].a_x1 = fa1;
        fluxes[_cell_right1].a_x0 = fa2;
        fluxes[_cell_right2].a_x0 = fa3;
        fluxes[_cell_left].u_x1 = fu1;
        fluxes[_cell_right1].u_x0 = fu2;
        fluxes[_cell_right2].u_x0 = fu3;
    }

private:
    struct NonlinearSystem : public INonlinearSystem6
    {
    public:
        NonlinearSystem(double mult1, double mult2, double mult3, double root4_area0_1, double root4_area0_2,
                        double root4_area0_3, double beta1, double beta2, double beta3, double root2_area0_1,
                        double root2_area0_2, double root2_area0_3)
            : _mult1(mult1), _root4_area0_1(root4_area0_1), _bet1(beta1), _root2_area0_1(root2_area0_1), _mult2(mult2),
              _root4_area0_2(root4_area0_2), _bet2(beta2), _root2_area0_2(root2_area0_2), _mult3(mult3),
              _root4_area0_3(root4_area0_3), _bet3(beta3), _root2_area0_3(root2_area0_3)
        {
            set_ww(0, 0, 0);
        }

        void set_ww(double w1, double w2, double w3)
        {
            _w1 = w1;
            _w2 = w2;
            _w3 = w3;
        };

        std::array<double, 6> f(double area1, double velo1, double area2, double velo2, double area3,
                                double velo3) const override
        {
            return {velo1 + _mult1 * (std::pow(area1, 0.25) - _root4_area0_1) - _w1,
                    velo2 - _mult2 * (std::pow(area2, 0.25) - _root4_area0_2) - _w2,
                    velo3 - _mult3 * (std::pow(area3, 0.25) - _root4_area0_3) - _w3,
                    area1 * velo1 - (area2 * velo2 + area3 * velo3),
                    velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo2 * velo2 / 2 -
                        _bet2 * (sqrt(area2) - _root2_area0_2),
                    velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo3 * velo3 / 2 -
                        _bet3 * (sqrt(area3) - _root2_area0_3)};
        };
        std::array<double, 36> jac(double area1, double velo1, double area2, double velo2, double area3,
                                   double velo3) const override
        {

            return {std::pow(area1, -0.75) * _mult1 / 4,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -std::pow(area2, -0.75) * _mult2 / 4,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -std::pow(area3, -0.75) * _mult3 / 4,
                    1.0,
                    velo1,
                    area1,
                    -velo2,
                    -area2,
                    -velo3,
                    -area3,
                    std::pow(area1, -0.5) * _bet1 / 2,
                    velo1,
                    -std::pow(area2, -0.5) * _bet2 / 2,
                    -velo2,
                    0.0,
                    0.0,
                    std::pow(area1, -0.5) * _bet1 / 2,
                    velo1,
                    0.0,
                    0.0,
                    -std::pow(area3, -0.5) * _bet3 / 2,
                    -velo3};
        };

    private:
        double _mult1, _root4_area0_1, _bet1, _root2_area0_1;
        double _mult2, _root4_area0_2, _bet2, _root2_area0_2;
        double _mult3, _root4_area0_3, _bet3, _root2_area0_3;
        double _w1, _w2, _w3;
    };
    const ProblemData& _data1;
    const ProblemData& _data2;
    const ProblemData& _data3;
    const size_t _node_left;
    const size_t _node_right1;
    const size_t _node_right2;
    const size_t _cell_left;
    const size_t _cell_right1;
    const size_t _cell_right2;

    NonlinearSystem _sys;

    double _area1_upw;
    double _area2_upw;
    double _area3_upw;
    double _velo1_upw;
    double _velo2_upw;
    double _velo3_upw;
};

class Junction3FluxCalculator : public IUpwindFluxCalculator
{
    Junction3FluxCalculator(const FemGrid& grid, const ProblemData& data1, const ProblemData& data2,
                            const ProblemData& data3, size_t cell_left1, size_t cell_left2, size_t cell_right)
        : _data1(data1), _data2(data2), _data3(data3), _node_left1(grid.tab_elem_nodes(cell_left1)[1]),
          _node_left2(grid.tab_elem_nodes(cell_left2)[1]), _node_right(grid.tab_elem_nodes(cell_right)[0]),
          _cell_left1(cell_left1), _cell_left2(cell_left2), _cell_right(cell_right),
          _sys(4 * sqrt(data1.beta / 2 / data1.rho), 4 * sqrt(data1.beta / 2 / data1.rho),
               4 * sqrt(data3.beta / 2 / data3.rho), data1.root4_a0, data2.root4_a0, data3.root4_a0, data1.beta,
               data2.beta, data3.beta, sqrt(data1.area0), sqrt(data2.area0), sqrt(data3.area0))
    {
        _area1_upw = _data1.area0;
        _area2_upw = _data2.area0;
        _area3_upw = _data3.area0;
        _velo1_upw = 0;
        _velo2_upw = 0;
        _velo3_upw = 0;
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double w1_upw = _data1.w1(area[_node_left1], velocity[_node_left1]);
        double w2_upw = _data2.w1(area[_node_left2], velocity[_node_left2]);
        double w3_upw = _data2.w2(area[_node_right], velocity[_node_right]);

        _sys.set_ww(w1_upw, w2_upw, w3_upw);
        solve_nonlinear_system(_sys, _area1_upw, _velo1_upw, _area2_upw, _velo2_upw, _area3_upw, _velo3_upw, 1e-12);

        double fa1 = _data1.flux_a(_area1_upw, _velo1_upw);
        double fu1 = _data1.flux_u(_area1_upw, _velo1_upw);
        double fa2 = _data2.flux_a(_area2_upw, _velo2_upw);
        double fu2 = _data2.flux_u(_area2_upw, _velo2_upw);
        double fa3 = _data2.flux_a(_area3_upw, _velo3_upw);
        double fu3 = _data2.flux_u(_area3_upw, _velo3_upw);

        fluxes[_cell_left1].a_x1 = fa1;
        fluxes[_cell_left2].a_x1 = fa2;
        fluxes[_cell_right].a_x0 = fa3;
        fluxes[_cell_left1].u_x1 = fu1;
        fluxes[_cell_left2].u_x1 = fu2;
        fluxes[_cell_right].u_x0 = fu3;
    }

private:
    struct NonlinearSystem : public INonlinearSystem6
    {
    public:
        NonlinearSystem(double mult1, double mult2, double mult3, double root4_area0_1, double root4_area0_2,
                        double root4_area0_3, double beta1, double beta2, double beta3, double root2_area0_1,
                        double root2_area0_2, double root2_area0_3)
            : _mult1(mult1), _root4_area0_1(root4_area0_1), _bet1(beta1), _root2_area0_1(root2_area0_1), _mult2(mult2),
              _root4_area0_2(root4_area0_2), _bet2(beta2), _root2_area0_2(root2_area0_2), _mult3(mult3),
              _root4_area0_3(root4_area0_3), _bet3(beta3), _root2_area0_3(root2_area0_3)
        {
            set_ww(0, 0, 0);
        }

        void set_ww(double w1, double w2, double w3)
        {
            _w1 = w1;
            _w2 = w2;
            _w3 = w3;
        };

        std::array<double, 6> f(double area1, double velo1, double area2, double velo2, double area3,
                                double velo3) const override
        {
            return {velo1 + _mult1 * (std::pow(area1, 0.25) - _root4_area0_1) - _w1,
                    velo2 + _mult2 * (std::pow(area2, 0.25) - _root4_area0_2) - _w2,
                    velo3 - _mult3 * (std::pow(area3, 0.25) - _root4_area0_3) - _w3,
                    area1 * velo1 + area2 * velo2 - area3 * velo3,
                    velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo2 * velo2 / 2 -
                        _bet2 * (sqrt(area2) - _root2_area0_2),
                    velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo3 * velo3 / 2 -
                        _bet3 * (sqrt(area3) - _root2_area0_3)};
        };
        std::array<double, 36> jac(double area1, double velo1, double area2, double velo2, double area3,
                                   double velo3) const override
        {

            return {std::pow(area1, -0.75) * _mult1 / 4,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    std::pow(area2, -0.75) * _mult2 / 4,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -std::pow(area3, -0.75) * _mult3 / 4,
                    1.0,
                    velo1,
                    area1,
                    velo2,
                    area2,
                    -velo3,
                    -area3,
                    std::pow(area1, -0.5) * _bet1 / 2,
                    velo1,
                    -std::pow(area2, -0.5) * _bet2 / 2,
                    -velo2,
                    0.0,
                    0.0,
                    std::pow(area1, -0.5) * _bet1 / 2,
                    velo1,
                    0.0,
                    0.0,
                    -std::pow(area3, -0.5) * _bet3 / 2,
                    -velo3};
        };

    private:
        double _mult1, _root4_area0_1, _bet1, _root2_area0_1;
        double _mult2, _root4_area0_2, _bet2, _root2_area0_2;
        double _mult3, _root4_area0_3, _bet3, _root2_area0_3;
        double _w1, _w2, _w3;
    };
    const ProblemData& _data1;
    const ProblemData& _data2;
    const ProblemData& _data3;
    const size_t _node_left1;
    const size_t _node_left2;
    const size_t _node_right;
    const size_t _cell_left1;
    const size_t _cell_left2;
    const size_t _cell_right;

    NonlinearSystem _sys;

    double _area1_upw;
    double _area2_upw;
    double _area3_upw;
    double _velo1_upw;
    double _velo2_upw;
    double _velo3_upw;
};
} 

TEST_CASE("Single vessel, inviscid", "[single-vessel-inviscid-explicit]")
{

    ProblemData data;
    double time = 0;
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {data.L};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.01, 1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1, nodes_coo);
    //FemGrid grid(data.L, 100, 1);
    double tau = grid.h() / 100;
    std::vector<ElementBoundaryFluxes> upwind_fluxes(grid.n_elements());

    NonstatGridSaver saver(grid1, nodes_coo, "bububu");
    saver.new_time_step(0);

    std::vector<double> velocity(grid.n_nodes(), 0.0);
    std::vector<double> area(grid.n_nodes(), data.area0);
    std::vector<double> pressure(grid.n_nodes(), 0.0);
    saver.save_vtk_point_data(velocity, "velocity");
    saver.save_vtk_point_data(area, "area");
    saver.save_vtk_point_data(pressure, "pressure");

    std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
    upwind_flux_calculator[0].reset(new InflowQFluxCalculator(
        grid, data, [&time]() { return time; }, 0));
    for (size_t i = 1; i < grid.n_points() - 1; ++i)
    {
        upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i - 1, i));
    }
    upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements() - 1));
    // prepare matrix solver
    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix tran = grid.transport_matrix();
    AmgcMatrixSolver slv;
    slv.set_matrix(mass);

    double t = 0.005;
    while (time < 0.2 - 1e-6)
    {
        time += tau;
        std::cout << "TIME=" << time;
        std::cout << "  Q=" << data.q_inflow(time) << std::endl;

        // nodewise fluxes
        std::vector<double> flux_a(grid.n_nodes());
        std::vector<double> flux_u(grid.n_nodes());
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            flux_a[i] = data.flux_a(area[i], velocity[i]);
            flux_u[i] = data.flux_u(area[i], velocity[i]);
        }

        // upwind fluxes
        for (auto c : upwind_flux_calculator)
        {
            c->compute(area, velocity, upwind_fluxes);
        }

        // assemble rhs
        // E
        std::vector<double> rhs_a = mass.mult_vec(area);
        std::vector<double> rhs_u = mass.mult_vec(velocity);
        // +tau*T
        std::vector<double> tran_a = tran.mult_vec(flux_a);
        std::vector<double> tran_u = tran.mult_vec(flux_u);
        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            rhs_a[i] -= tau * tran_a[i];
            rhs_u[i] -= tau * tran_u[i];
        }
        // + coupling
        for (size_t ielem = 0; ielem < grid.n_elements(); ++ielem)
        {
            size_t node0 = grid.tab_elem_nodes(ielem)[0];
            size_t node1 = grid.tab_elem_nodes(ielem)[1];
            rhs_a[node0] += tau * upwind_fluxes[ielem].a_x0;
            rhs_a[node1] -= tau * upwind_fluxes[ielem].a_x1;
            rhs_u[node0] += tau * upwind_fluxes[ielem].u_x0;
            rhs_u[node1] -= tau * upwind_fluxes[ielem].u_x1;
        }
        // solve
        slv.solve(rhs_a, area);
        slv.solve(rhs_u, velocity);

        for (size_t i = 0; i < grid.n_nodes(); ++i)
        {
            pressure[i] = data.pressure(area[i]);
        }
        if (t - 1e-5< time )
        {
            saver.new_time_step(time);
            saver.save_vtk_point_data(velocity, "velocity");
            saver.save_vtk_point_data(area, "area");
            saver.save_vtk_point_data(pressure, "pressure");
            t += 0.005;
        }
    }

    double maxp = *std::max_element(pressure.begin(), pressure.end());
    CHECK(maxp == Approx(24.3958).margin(1e-3));
}
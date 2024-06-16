#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP
#include "fem_grid.hpp"
#include "tests/nonlinsolver.hpp"
#include <functional>

namespace bflow
{
struct ElementBoundaryFluxes
{
    double upwind_area_x0 = 0;
    double upwind_area_x1 = 0;
    double upwind_velo_x0 = 0;
    double upwind_velo_x1 = 0;
};

struct ProblemData
{
    static constexpr double pi = 3.1415926;

    ProblemData();
    void recompute();
    void report() const;
    // geometry parameters
    double area0 = pi * 1e-4;
    // fluid parameters
    double rho = 1050;
    double mu = 0;
    double profile_order = 9;
    // vessel tissue parameters
    double h = 1.5e-3;
    double E = 4e5;
    // p(a)
    double pressure(double area) const;
    // fluxes
    double flux_a(double area, double vel) const;
    double flux_u(double area, double vel) const;
    // characteristic variables
    double w1(double area, double vel) const;
    double w2(double area, double vel) const;
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
    InflowQFluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> getq, size_t cell,
                          double eps = 1e-12);
    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

private:
    struct NonlinearSystem : public INonlinearSystem2
    {
    public:
        NonlinearSystem(double mult, double root4_area0);
        void set_qw(double q, double w2);
        std::array<double, 2> f(double area, double velo) const override;
        std::array<double, 4> jac(double area, double velo) const override;

    private:
        const double _mult, _root4_area0;
        double _q, _w2;
    };

    const ProblemData& _data;
    std::function<double()> _getq;
    const size_t _cell_right;
    const size_t _node_right;
    NonlinearSystem _sys;
    const double _eps;

    double _a;
};
class InflowPFluxCalculator_NonReflecting : public IUpwindFluxCalculator
{
public:
    InflowPFluxCalculator_NonReflecting(const FemGrid& grid, const ProblemData& data, std::function<double()> getp,
                                        size_t cell)
        : _data(data), _getp(getp), _cell_right(cell), _node_right(grid.tab_elem_nodes(cell)[0])
    {
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

private:
    const ProblemData& _data;
    std::function<double()> _getp;
    const size_t _cell_right;
    const size_t _node_right;
};

class InflowPFluxCalculator : public IUpwindFluxCalculator
{
public:
    InflowPFluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> getp, size_t cell,
                          double eps = 1e-12)
        : _data(data), _getp(getp), _cell_right(cell), _node_right(grid.tab_elem_nodes(cell)[0]),
          _sys(data.amult, data.area0, data.beta), _eps(eps)
    {
        _a = data.area0;
        _v = 0;
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

private:
    struct NonlinearSystem : public INonlinearSystem2
    {
    public:
        NonlinearSystem(double mult, double area0, double beta)
            : _mult(mult), _area0(area0), _beta(beta), _root2_area0(std::pow(area0, 0.5)),
              _root4_area0(std::pow(area0, 0.25))
        {
            set_pw(0, 0);
        }

        void set_pw(double p, double w2)
        {
            _p = p;
            _w2 = w2;
        };

        std::array<double, 2> f(double area, double velo) const override
        {
            return {_beta * (sqrt(area) - _root2_area0) - _p, velo - _mult * (sqrt(sqrt(area)) - _root4_area0) - _w2};
        };
        std::array<double, 4> jac(double area, double velo) const override
        {
            return {0.5 * _beta / sqrt(area), 0, -0.25 * _mult * pow(area, -0.75), 1};
        };

    private:
        const double _mult, _area0, _beta, _root2_area0, _root4_area0;
        double _p, _w2;
    };

    const ProblemData& _data;
    std::function<double()> _getp;
    const size_t _cell_right;
    const size_t _node_right;
    NonlinearSystem _sys;
    const double _eps;
    double _a, _v;
};

class InflowW1FluxCalculator : public IUpwindFluxCalculator
{
public:
    InflowW1FluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> getw1, size_t cell)
        : _data(data), _getw1(getw1), _cell_right(cell), _node_right(grid.tab_elem_nodes(cell)[0])
    {
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

private:
    const ProblemData& _data;
    std::function<double()> _getw1;
    const size_t _cell_right;
    const size_t _node_right;
};
class InflowAUFluxCalculator : public IUpwindFluxCalculator
{
public:
    InflowAUFluxCalculator(const FemGrid& grid, const ProblemData& data,
                           std::function<std::pair<double, double>()> getau, size_t cell)
        : _data(data), _cell_right(cell), _getau(getau)
    {
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double a, u;
        std::tie(a, u) = _getau();
        // fluxes[_cell_right].upwind_area_x0 = a;
        // fluxes[_cell_right].upwind_velo_x0 = u;

        double w1_upw = _data.w1(a, u);
        size_t _node_right = 0;
        double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);
        // double w2_upw = 0;
        double a1 = (w1_upw - w2_upw) / 2 / _data.amult + _data.root4_a0;
        double area_upw = a1 * a1 * a1 * a1;
        double velo_upw = (w1_upw + w2_upw) / 2.0;

        fluxes[_cell_right].upwind_area_x0 = area_upw;
        fluxes[_cell_right].upwind_velo_x0 = velo_upw;
    }

private:
    const ProblemData& _data;
    const size_t _cell_right;
    std::function<std::pair<double, double>()> _getau;
};

class OutflowFluxCalculator : public IUpwindFluxCalculator
{
public:
    OutflowFluxCalculator(const FemGrid& grid, const ProblemData& data, size_t cell)
        : _data(data), _cell_left(cell), _node_left(grid.tab_elem_nodes(cell)[1])
    {
    }
    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

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
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

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
    MergingFluxCalculator(const FemGrid& grid, const ProblemData& data_left, const ProblemData& data_right,
                          size_t cell_left, size_t cell_right, double eps = 1e-12);
    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

private:
    struct NonlinearSystem : public INonlinearSystem4
    {
    public:
        NonlinearSystem(const ProblemData& data1, const ProblemData& data2);
        void set_ww(double w1, double w2);
        std::array<double, 4> f(double area1, double velo1, double area2, double velo2) const override;
        std::array<double, 16> jac(double area1, double velo1, double area2, double velo2) const override;

    private:
        const ProblemData& _data1;
        const ProblemData& _data2;
        double _w1, _w2;
    };
    const ProblemData& _data_left;
    const ProblemData& _data_right;
    const size_t _node_left, _node_right;
    const size_t _cell_left, _cell_right;
    NonlinearSystem _sys;
    const double _eps;

    double _a1, _a2, _u1, _u2;
};

class Bifurcation3FluxCalculator : public IUpwindFluxCalculator
{
public:
    Bifurcation3FluxCalculator(const FemGrid& grid, const ProblemData& data_left, const ProblemData& data_right1,
                               const ProblemData& data_right2, size_t elem_left, size_t elem_right1, size_t elem_right2,
                               double eps = 1e-12)
        : _elem1(elem_left), _elem2(elem_right1), _elem3(elem_right2), _node1(grid.tab_elem_nodes(elem_left)[1]),
          _node2(grid.tab_elem_nodes(elem_right1)[0]), _node3(grid.tab_elem_nodes(elem_right2)[0]), _data1(data_left),
          _data2(data_right1), _data3(data_right2), _eps(eps), _sys(_data1, _data2, _data3), _a1(_data1.area0),
          _a2(_data2.area0), _a3(_data3.area0), _u1(0), _u2(0), _u3(0)
    {
    }

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override
    {
        double w1 = _data1.w1(area[_node1], velocity[_node1]);
        double w2 = _data2.w2(area[_node2], velocity[_node2]);
        double w3 = _data3.w2(area[_node3], velocity[_node3]);

        double area1 = _a1;
        double velo1 = _u1;
        double area2 = _a2;
        double velo2 = _u2;
        double area3 = _a3;
        double velo3 = _u3;

        _sys.set_www(w1, w2, w3);
        solve_nonlinear_system(_sys, area1, velo1, area2, velo2, area3, velo3, _eps);

        fluxes[_elem1].upwind_area_x1 = area1;
        fluxes[_elem1].upwind_velo_x1 = velo1;
        fluxes[_elem2].upwind_area_x0 = area2;
        fluxes[_elem2].upwind_velo_x0 = velo2;
        fluxes[_elem3].upwind_area_x0 = area3;
        fluxes[_elem3].upwind_velo_x0 = velo3;

        _a1 = area1;
        _u1 = velo1;
        _a2 = area2;
        _u2 = velo2;
        _a3 = area3;
        _u3 = velo3;
    }

public:
    struct NonlinearSystem : public INonlinearSystem6
    {

        NonlinearSystem(const ProblemData& data1, const ProblemData& data2, const ProblemData& data3)
            : _data1(data1), _data2(data2), _data3(data3)
        {
        }

        std::array<double, 6> f(double area1, double velo1, double area2, double velo2, double area3,
                                double velo3) const override
        {
            return {_data1.flux_a(area1, velo1) - _data2.flux_a(area2, velo2) - _data3.flux_a(area3, velo3),
                    _data1.flux_u(area1, velo1) - _data2.flux_u(area2, velo2),
                    _data1.flux_u(area1, velo1) - _data3.flux_u(area3, velo3),
                    _data1.w1(area1, velo1) - _w1,
                    _data2.w2(area2, velo2) - _w2,
                    _data3.w2(area3, velo3) - _w3};
        }
        std::array<double, 36> jac(double area1, double velo1, double area2, double velo2, double area3,
                                   double velo3) const override
        {
            return {velo1,
                    area1,
                    -velo2,
                    -area2,
                    -velo3,
                    -area3,
                    0.5 * _data1.beta / _data1.rho / sqrt(area1),
                    velo1,
                    -0.5 * _data2.beta / _data2.rho / sqrt(area2),
                    -velo2,
                    0,
                    0,
                    0.5 * _data1.beta / _data1.rho / sqrt(area1),
                    velo1,
                    0,
                    0,
                    -0.5 * _data3.beta / _data3.rho / sqrt(area3),
                    -velo3,
                    sqrt(_data1.beta / 2 / _data1.rho) * std::pow(area1, -0.75),
                    1,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    -sqrt(_data2.beta / 2 / _data2.rho) * std::pow(area2, -0.75),
                    1,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    -sqrt(_data3.beta / 2 / _data3.rho) * std::pow(area3, -0.75),
                    1};
        }

        void set_www(double w1, double w2, double w3)
        {
            _w1 = w1;
            _w2 = w2;
            _w3 = w3;
        }

    private:
        const ProblemData& _data1;
        const ProblemData& _data2;
        const ProblemData& _data3;

        double _w1 = 0, _w2 = 0, _w3 = 0;
    };

    const size_t _elem1, _elem2, _elem3;
    const size_t _node1, _node2, _node3;
    const ProblemData &_data1, _data2, _data3;
    const double _eps;
    NonlinearSystem _sys;

    double _a1, _a2, _a3;
    double _u1, _u2, _u3;
};

class Junction3FluxCalculator : public IUpwindFluxCalculator
{
    Junction3FluxCalculator(const FemGrid& grid, const ProblemData& data1, const ProblemData& data2,
                            const ProblemData& data3, size_t cell_left1, size_t cell_left2, size_t cell_right);

    void compute(const std::vector<double>& area, const std::vector<double>& velocity,
                 std::vector<ElementBoundaryFluxes>& fluxes) override;

private:
    struct NonlinearSystem : public INonlinearSystem6
    {
    public:
        NonlinearSystem(double mult1, double mult2, double mult3, double root4_area0_1, double root4_area0_2,
                        double root4_area0_3, double beta1, double beta2, double beta3, double root2_area0_1,
                        double root2_area0_2, double root2_area0_3);
        void set_ww(double w1, double w2, double w3);
        std::array<double, 6> f(double area1, double velo1, double area2, double velo2, double area3,
                                double velo3) const override;
        std::array<double, 36> jac(double area1, double velo1, double area2, double velo2, double area3,
                                   double velo3) const override;

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
} // namespace bflow
#endif

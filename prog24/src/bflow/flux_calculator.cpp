#include "bflow/flux_calculator.hpp"
using namespace bflow;

ProblemData::ProblemData()
{
    recompute();
}

void ProblemData::recompute()
{
    beta = 4.0 / 3.0 * sqrt(pi) * h * E / area0;
    amult = 4 * sqrt(beta / 2 / rho);
    root4_a0 = sqrt(sqrt(area0));
    visc_coef = -2 * (profile_order + 2) * mu * pi / profile_order / rho;
}
// inflow conditions
double ProblemData::q_inflow(double t) const
{
    return 1e-6 * exp(-1e4 * (t - 0.05) * (t - 0.05));
};

// p(a)
double ProblemData::pressure(double area) const
{
    return 0 + beta * (std::sqrt(area) - std::sqrt(area0));
};

// fluxes
double ProblemData::flux_a(double area, double vel) const
{
    return area * vel;
};
double ProblemData::flux_u(double area, double vel) const
{
    return vel * vel / 2 + pressure(area) / rho;
}

// characteristic variables
double ProblemData::w1(double area, double vel) const
{
    return vel + amult * (sqrt(sqrt(area)) - root4_a0);
}
double ProblemData::w2(double area, double vel) const
{
    return vel - amult * (sqrt(sqrt(area)) - root4_a0);
}

InflowQFluxCalculator::InflowQFluxCalculator(const FemGrid& grid, const ProblemData& data,
                                             std::function<double()> get_time, size_t cell)
    : _data(data), _get_time(get_time), _cell_right(cell), _node_right(grid.tab_elem_nodes(cell)[0]),
      _sys(data.amult, data.root4_a0)
{
    _a = data.area0;
}

void InflowQFluxCalculator::compute(const std::vector<double>& area, const std::vector<double>& velocity,
                                    std::vector<ElementBoundaryFluxes>& fluxes)
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

    fluxes[_cell_right].upwind_area_x0 = area_upw;
    fluxes[_cell_right].upwind_velo_x0 = velo_upw;

    _a = area_upw;
}

InflowQFluxCalculator::NonlinearSystem::NonlinearSystem(double mult, double root4_area0)
    : _mult(mult), _root4_area0(root4_area0)
{
    set_qw(0, 0);
}

void InflowQFluxCalculator::NonlinearSystem::set_qw(double q, double w2)
{
    _q = q;
    _w2 = w2;
}

std::array<double, 2> InflowQFluxCalculator::NonlinearSystem::f(double area, double velo) const
{
    return {area * velo - _q, velo - _mult * (sqrt(sqrt(area)) - _root4_area0) - _w2};
}

std::array<double, 4> InflowQFluxCalculator::NonlinearSystem::jac(double area, double velo) const
{
    return {velo, area, -0.25 * _mult * pow(area, -0.75), 1};
}

void OutflowFluxCalculator::compute(const std::vector<double>& area, const std::vector<double>& velocity,
                                    std::vector<ElementBoundaryFluxes>& fluxes)
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

    fluxes[_cell_left].upwind_area_x1 = area_upw;
    fluxes[_cell_left].upwind_velo_x1 = velo_upw;
}

void InternalFluxCalculator::compute(const std::vector<double>& area, const std::vector<double>& velocity,
                                     std::vector<ElementBoundaryFluxes>& fluxes)
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

    
    fluxes[_cell_right].upwind_area_x0 = area_upw;
    fluxes[_cell_right].upwind_velo_x0 = velo_upw;
    fluxes[_cell_left].upwind_area_x1 = area_upw;
    fluxes[_cell_left].upwind_velo_x1 = velo_upw;
}

MergingFluxCalculator::MergingFluxCalculator(const FemGrid& grid, const ProblemData& datal, const ProblemData& datar,
                                             size_t cell_left, size_t cell_right)
    : _data1(datal), _data2(datar), _node_left(grid.tab_elem_nodes(cell_left)[1]),
      _node_right(grid.tab_elem_nodes(cell_right)[0]), _cell_left(cell_left), _cell_right(cell_right),
      _sys(4 * sqrt(datal.beta / 2 / datal.rho), 4 * sqrt(datar.beta / 2 / datar.rho), datal.root4_a0, datar.root4_a0,
           datal.beta, datar.beta, sqrt(datal.area0), sqrt(datar.area0))
{
    _area1_upw = _data1.area0;
    _area2_upw = _data2.area0;
    _velo1_upw = 0;
    _velo2_upw = 0;
}

void MergingFluxCalculator::compute(const std::vector<double>& area, const std::vector<double>& velocity,
                                    std::vector<ElementBoundaryFluxes>& fluxes)
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

MergingFluxCalculator::NonlinearSystem::NonlinearSystem(double mult1, double mult2, double root4_area0_1,
                                                        double root4_area0_2, double beta1, double beta2,
                                                        double root2_area0_1, double root2_area0_2)
    : _mult1(mult1), _root4_area0_1(root4_area0_1), _bet1(beta1), _root2_area0_1(root2_area0_1), _mult2(mult2),
      _root4_area0_2(root4_area0_2), _bet2(beta2), _root2_area0_2(root2_area0_2)
{
    set_ww(0, 0);
}

void MergingFluxCalculator::NonlinearSystem::set_ww(double w1, double w2)
{
    _w1 = w1;
    _w2 = w2;
};

std::array<double, 4> MergingFluxCalculator::NonlinearSystem::f(double area1, double velo1, double area2,
                                                                double velo2) const
{
    // clang-format off
    return {velo1 + _mult1 * (std::pow(area1, 0.25) - _root4_area0_1) - _w1,
            velo2 - _mult2 * (std::pow(area2, 0.25) - _root4_area0_2) - _w2,
            area1 * velo1 - area2 * velo2,
            velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo2 * velo2 / 2 - _bet2 * (sqrt(area2) - _root2_area0_2)};
    // clang-format on
}
std::array<double, 16> MergingFluxCalculator::NonlinearSystem::jac(double area1, double velo1, double area2,
                                                                   double velo2) const
{
    // clang-format off
    return {std::pow(area1, -0.75) * _mult1 / 4, 1.0, 0.0, 0.0,
            0.0, 0.0, -std::pow(area2, -0.75) * _mult2 / 4, 1.0,
            velo1, area1, -velo2, -area2,
            std::pow(area1, -0.5) * _bet1 / 2, velo1, -std::pow(area2, -0.5) * _bet2 / 2, -velo2};
    // clang-format on
}

Bifurcation3FluxCalculator::Bifurcation3FluxCalculator(const FemGrid& grid, const ProblemData& data1,
                                                       const ProblemData& data2, const ProblemData& data3,
                                                       size_t cell_left, size_t cell_right1, size_t cell_right2)
    : _data1(data1), _data2(data2), _data3(data3), _node_left(grid.tab_elem_nodes(cell_left)[1]),
      _node_right1(grid.tab_elem_nodes(cell_right1)[0]), _node_right2(grid.tab_elem_nodes(cell_right2)[0]),
      _cell_left(cell_left), _cell_right1(cell_right1), _cell_right2(cell_right2),
      _sys(4 * sqrt(data1.beta / 2 / data1.rho), 4 * sqrt(data1.beta / 2 / data1.rho),
           4 * sqrt(data3.beta / 2 / data3.rho), data1.root4_a0, data2.root4_a0, data3.root4_a0, data1.beta, data2.beta,
           data3.beta, sqrt(data1.area0), sqrt(data2.area0), sqrt(data3.area0))

{
    _area1_upw = _data1.area0;
    _area2_upw = _data2.area0;
    _area3_upw = _data3.area0;
    _velo1_upw = 0;
    _velo2_upw = 0;
    _velo3_upw = 0;
}

void Bifurcation3FluxCalculator::compute(const std::vector<double>& area, const std::vector<double>& velocity,
                                         std::vector<ElementBoundaryFluxes>& fluxes)
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

Bifurcation3FluxCalculator::NonlinearSystem::NonlinearSystem(double mult1, double mult2, double mult3,
                                                             double root4_area0_1, double root4_area0_2,
                                                             double root4_area0_3, double beta1, double beta2,
                                                             double beta3, double root2_area0_1, double root2_area0_2,
                                                             double root2_area0_3)
    : _mult1(mult1), _root4_area0_1(root4_area0_1), _bet1(beta1), _root2_area0_1(root2_area0_1), _mult2(mult2),
      _root4_area0_2(root4_area0_2), _bet2(beta2), _root2_area0_2(root2_area0_2), _mult3(mult3),
      _root4_area0_3(root4_area0_3), _bet3(beta3), _root2_area0_3(root2_area0_3)
{
    set_ww(0, 0, 0);
}

void Bifurcation3FluxCalculator::NonlinearSystem::set_ww(double w1, double w2, double w3)
{
    _w1 = w1;
    _w2 = w2;
    _w3 = w3;
};

std::array<double, 6> Bifurcation3FluxCalculator::NonlinearSystem::f(double area1, double velo1, double area2,
                                                                     double velo2, double area3, double velo3) const
{

    // clang-format off
    return {velo1 + _mult1 * (std::pow(area1, 0.25) - _root4_area0_1) - _w1,
            velo2 - _mult2 * (std::pow(area2, 0.25) - _root4_area0_2) - _w2,
            velo3 - _mult3 * (std::pow(area3, 0.25) - _root4_area0_3) - _w3,
            area1 * velo1 - (area2 * velo2 + area3 * velo3),
            velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo2 * velo2 / 2 - _bet2 * (sqrt(area2) - _root2_area0_2),
            velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo3 * velo3 / 2 - _bet3 * (sqrt(area3) - _root2_area0_3)};
    // clang-format on
};
std::array<double, 36> Bifurcation3FluxCalculator::NonlinearSystem::jac(double area1, double velo1, double area2,
                                                                        double velo2, double area3, double velo3) const
{
    // clang-format off
    return {std::pow(area1, -0.75) * _mult1 / 4, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, -std::pow(area2, -0.75) * _mult2 / 4, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, -std::pow(area3, -0.75) * _mult3 / 4, 1.0,
            velo1, area1, -velo2, -area2, -velo3, -area3,
            std::pow(area1, -0.5) * _bet1 / 2, velo1, -std::pow(area2, -0.5) * _bet2 / 2, -velo2, 0.0, 0.0,
            std::pow(area1, -0.5) * _bet1 / 2,velo1, 0.0, 0.0, -std::pow(area3, -0.5) * _bet3 / 2, -velo3};
    // clang-format on
};

Junction3FluxCalculator::Junction3FluxCalculator(const FemGrid& grid, const ProblemData& data1,
                                                 const ProblemData& data2, const ProblemData& data3, size_t cell_left1,
                                                 size_t cell_left2, size_t cell_right)
    : _data1(data1), _data2(data2), _data3(data3), _node_left1(grid.tab_elem_nodes(cell_left1)[1]),
      _node_left2(grid.tab_elem_nodes(cell_left2)[1]), _node_right(grid.tab_elem_nodes(cell_right)[0]),
      _cell_left1(cell_left1), _cell_left2(cell_left2), _cell_right(cell_right),
      _sys(4 * sqrt(data1.beta / 2 / data1.rho), 4 * sqrt(data1.beta / 2 / data1.rho),
           4 * sqrt(data3.beta / 2 / data3.rho), data1.root4_a0, data2.root4_a0, data3.root4_a0, data1.beta, data2.beta,
           data3.beta, sqrt(data1.area0), sqrt(data2.area0), sqrt(data3.area0))
{
    _area1_upw = _data1.area0;
    _area2_upw = _data2.area0;
    _area3_upw = _data3.area0;
    _velo1_upw = 0;
    _velo2_upw = 0;
    _velo3_upw = 0;
}

void Junction3FluxCalculator::compute(const std::vector<double>& area, const std::vector<double>& velocity,
                                      std::vector<ElementBoundaryFluxes>& fluxes)
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

Junction3FluxCalculator::NonlinearSystem::NonlinearSystem(double mult1, double mult2, double mult3,
                                                          double root4_area0_1, double root4_area0_2,
                                                          double root4_area0_3, double beta1, double beta2,
                                                          double beta3, double root2_area0_1, double root2_area0_2,
                                                          double root2_area0_3)
    : _mult1(mult1), _root4_area0_1(root4_area0_1), _bet1(beta1), _root2_area0_1(root2_area0_1), _mult2(mult2),
      _root4_area0_2(root4_area0_2), _bet2(beta2), _root2_area0_2(root2_area0_2), _mult3(mult3),
      _root4_area0_3(root4_area0_3), _bet3(beta3), _root2_area0_3(root2_area0_3)
{
    set_ww(0, 0, 0);
}

void Junction3FluxCalculator::NonlinearSystem::set_ww(double w1, double w2, double w3)
{
    _w1 = w1;
    _w2 = w2;
    _w3 = w3;
};

std::array<double, 6> Junction3FluxCalculator::NonlinearSystem::f(double area1, double velo1, double area2,
                                                                  double velo2, double area3, double velo3) const
{
    // clang-format off
    return {velo1 + _mult1 * (std::pow(area1, 0.25) - _root4_area0_1) - _w1,
            velo2 + _mult2 * (std::pow(area2, 0.25) - _root4_area0_2) - _w2,
            velo3 - _mult3 * (std::pow(area3, 0.25) - _root4_area0_3) - _w3,
            area1 * velo1 + area2 * velo2 - area3 * velo3,
            velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo2 * velo2 / 2 - _bet2 * (sqrt(area2) - _root2_area0_2),
            velo1 * velo1 / 2 + _bet1 * (sqrt(area1) - _root2_area0_1) - velo3 * velo3 / 2 - _bet3 * (sqrt(area3) - _root2_area0_3)};
    // clang-format on
};
std::array<double, 36> Junction3FluxCalculator::NonlinearSystem::jac(double area1, double velo1, double area2,
                                                                     double velo2, double area3, double velo3) const
{
    // clang-format off
    return {std::pow(area1, -0.75) * _mult1 / 4, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, std::pow(area2, -0.75) * _mult2 / 4, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, -std::pow(area3, -0.75) * _mult3 / 4, 1.0,
            velo1, area1, velo2, area2, -velo3, -area3,
            std::pow(area1, -0.5) * _bet1 / 2, velo1, -std::pow(area2, -0.5) * _bet2 / 2, -velo2, 0.0, 0.0,
            std::pow(area1, -0.5) * _bet1 / 2, velo1, 0.0, 0.0, -std::pow(area3, -0.5) * _bet3 / 2, -velo3};
    // clang-format on
};
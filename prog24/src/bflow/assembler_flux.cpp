#include "bflow/assembler_flux.hpp"
using namespace bflow;

AssemblerFlux::AssemblerFlux(const FemGrid& grid, const ProblemData& data) : _grid(grid), _data(data)
{
    _mass = _grid.mass_matrix();
    _tran = _grid.transport_matrix();

    _upwind_flux_calculator.resize(grid.n_points());
    _upwind_flux_calculator[0].reset(new InflowQFluxCalculator(
        grid, data, [&]() { return _time; }, 0));
    for (size_t i = 1; i < grid.n_points() - 1; ++i)
    {
        _upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i - 1, i));
    }
    _upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements() - 1));

    _upwind_fluxes.resize(grid.n_elements());
    _flux_a.resize(_grid.n_nodes());
    _flux_u.resize(_grid.n_nodes());
    _visc.resize(_grid.n_nodes());

    _load_vector = _mass.mult_vec(std::vector<double>(_grid.n_nodes(), 1.0));
}

const CsrMatrix& AssemblerFlux::mass() const
{
    return _mass;
}

const CsrMatrix& AssemblerFlux::block_transport() const
{
    return _tran;
}

void AssemblerFlux::actualize_fluxes(double time, const std::vector<double>& area, const std::vector<double>& velo)
{
    _time = time;
    _u = velo;
    _area = area;

    // compute upwind fluxes
    for (auto c : _upwind_flux_calculator)
    {
        c->compute(area, velo, _upwind_fluxes);
    }

    // nodewise fluxes
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        _flux_a[i] = _data.flux_a(area[i], velo[i]);
        _flux_u[i] = _data.flux_u(area[i], velo[i]);
    }

    // viscous term
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        _visc[i] = _data.visc_coef * velo[i] / area[i];
    }
}

std::vector<double> AssemblerFlux::dfa_dx() const
{
    // nodewise flux
    std::vector<double> ret = _tran.mult_vec(_flux_a);
    for (auto& v : ret)
        v *= -1;

    // coupling
    for (size_t ielem = 0; ielem < _grid.n_elements(); ++ielem)
    {
        size_t node0 = _grid.tab_elem_nodes(ielem)[0];
        size_t node1 = _grid.tab_elem_nodes(ielem)[1];
        ret[node0] -= _upwind_fluxes[ielem].a_x0;
        ret[node1] += _upwind_fluxes[ielem].a_x1;
    }

    return ret;
}
std::vector<double> AssemblerFlux::dfu_dx() const
{
    // nodewise flux
    std::vector<double> ret = _tran.mult_vec(_flux_u);
    for (auto& v : ret)
        v *= -1;

    // coupling
    for (size_t ielem = 0; ielem < _grid.n_elements(); ++ielem)
    {
        size_t node0 = _grid.tab_elem_nodes(ielem)[0];
        size_t node1 = _grid.tab_elem_nodes(ielem)[1];
        ret[node0] -= _upwind_fluxes[ielem].u_x0;
        ret[node1] += _upwind_fluxes[ielem].u_x1;
    }

    return ret;
}

std::vector<double> AssemblerFlux::viscous_term() const
{
    return _mass.mult_vec(_visc);
}

CsrMatrix AssemblerFlux::block_u_transport() const
{
    CsrMatrix ret(_tran);
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        for (size_t a = ret.addr()[i]; a < ret.addr()[i + 1]; ++a)
        {
            size_t col = ret.cols()[a];
            ret.vals()[a] *= _u[col];
        }
    }
    return ret;
}

std::vector<double> AssemblerFlux::coupling_flux_ua() const
{
    std::vector<double> ret(_grid.n_nodes(), 0.0);
    for (size_t ielem = 0; ielem < _grid.n_elements(); ++ielem)
    {
        size_t node0 = _grid.tab_elem_nodes(ielem)[0];
        size_t node1 = _grid.tab_elem_nodes(ielem)[1];
        ret[node0] -= _upwind_fluxes[ielem].upwind_area_x0 * _upwind_fluxes[ielem].upwind_velo_x0;
        ret[node1] += _upwind_fluxes[ielem].upwind_area_x1 * _upwind_fluxes[ielem].upwind_velo_x1;
    }
    return ret;
}

std::vector<double> AssemblerFlux::coupling_flux_u2() const
{
    std::vector<double> ret(_grid.n_nodes(), 0.0);
    for (size_t ielem = 0; ielem < _grid.n_elements(); ++ielem)
    {
        size_t node0 = _grid.tab_elem_nodes(ielem)[0];
        size_t node1 = _grid.tab_elem_nodes(ielem)[1];
        ret[node0] -= _upwind_fluxes[ielem].upwind_velo_x0 * _upwind_fluxes[ielem].upwind_velo_x0;
        ret[node1] += _upwind_fluxes[ielem].upwind_velo_x1 * _upwind_fluxes[ielem].upwind_velo_x1;
    }
    return ret;
}

std::vector<double> AssemblerFlux::coupling_flux_p() const
{
    std::vector<double> ret(_grid.n_nodes(), 0.0);
    for (size_t ielem = 0; ielem < _grid.n_elements(); ++ielem)
    {
        size_t node0 = _grid.tab_elem_nodes(ielem)[0];
        size_t node1 = _grid.tab_elem_nodes(ielem)[1];
        ret[node0] -= _data.pressure(_upwind_fluxes[ielem].upwind_area_x0);
        ret[node1] += _data.pressure(_upwind_fluxes[ielem].upwind_area_x1);
    }
    return ret;
}

double AssemblerFlux::compute_residual(const CsrMatrix& lhs, const std::vector<double>& rhs,
                                       const std::vector<double>& u) const
{
    std::vector<double> r = lhs.mult_vec(u);
    for (size_t i = 0; i < r.size(); ++i)
    {
        r[i] = rhs[i] - r[i];
    }
    return vector_norm2(r);
}

double AssemblerFlux::vector_norm2(const std::vector<double>& v) const
{
    double s = 0;
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        s += v[i] * v[i] * _load_vector[i];
    }
    return std::sqrt(s / _grid.full_length());
}
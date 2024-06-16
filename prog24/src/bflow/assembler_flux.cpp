#include "assembler_flux.hpp"
using namespace bflow;

AssemblerFlux::AssemblerFlux(const FemGrid& grid, const std::vector<const ProblemData*>& data,
                             const std::vector<std::shared_ptr<IUpwindFluxCalculator>>& upwind_calculators)
    : _grid(grid), _data(data), _upwind_flux_calculator(upwind_calculators)
{
    _data_by_node = std::vector<const ProblemData*>(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        size_t icell = grid.tab_node_elem(i);
        _data_by_node[i] = _data[icell];
    }

    _mass = _grid.mass_matrix();
    _tran = _grid.block_transport_matrix();

    _upwind_fluxes.resize(grid.n_elements());
    _flux_a.resize(_grid.n_nodes());
    _flux_u.resize(_grid.n_nodes());

    _load_vector = _mass.mult_vec(std::vector<double>(_grid.n_nodes(), 1.0));
}
AssemblerFlux::AssemblerFlux(const FemGrid& grid, const ProblemData* data,
                             const std::vector<std::shared_ptr<IUpwindFluxCalculator>>& upwind_calculators)
    : AssemblerFlux(grid, std::vector<const ProblemData*>(grid.n_elements(), data), upwind_calculators)
{
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
        _flux_a[i] = _data_by_node[i]->flux_a(area[i], velo[i]);
        _flux_u[i] = _data_by_node[i]->flux_u(area[i], velo[i]);
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
        const auto& uf = _upwind_fluxes[ielem];
        ret[node0] -= uf.upwind_area_x0 * uf.upwind_velo_x0;
        ret[node1] += uf.upwind_area_x1 * uf.upwind_velo_x1;
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
        const auto& uf = _upwind_fluxes[ielem];

        ret[node0] -=
            0.5 * uf.upwind_velo_x0 * uf.upwind_velo_x0 + _data[ielem]->pressure(uf.upwind_area_x0) / _data[ielem]->rho;
        ret[node1] +=
            0.5 * uf.upwind_velo_x1 * uf.upwind_velo_x1 + _data[ielem]->pressure(uf.upwind_area_x1) / _data[ielem]->rho;
    }
    return ret;
}

CsrMatrix AssemblerFlux::viscous_matrix() const
{
    CsrMatrix ret(_mass);
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        for (size_t a = ret.addr()[i]; a < ret.addr()[i + 1]; ++a)
        {
            size_t col = ret.cols()[a];
            ret.vals()[a] *= (_data_by_node[i]->visc_coef / _area[col]);
        }
    }
    return ret;
}

CsrMatrix AssemblerFlux::block_u_transport() const
{
    // return _grid.block_u_transport_matrix(_u);
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
        ret[node0] -= _data[ielem]->pressure(_upwind_fluxes[ielem].upwind_area_x0);
        ret[node1] += _data[ielem]->pressure(_upwind_fluxes[ielem].upwind_area_x1);
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

std::vector<double> AssemblerFlux::pressure(const std::vector<double>& area) const
{
    const double* a = (area.size() == 0) ? _area.data() : area.data();
    std::vector<double> ret(_grid.n_nodes());
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        ret[i] = _data_by_node[i]->pressure(a[i]);
    }
    return ret;
};

std::vector<double> AssemblerFlux::w1() const
{
    const double* a = _area.data();
    const double* u = _u.data();
    std::vector<double> ret(_grid.n_nodes());
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        ret[i] = _data_by_node[i]->w1(a[i], u[i]);
    }
    return ret;
};

std::vector<double> AssemblerFlux::w2() const
{
    const double* a = _area.data();
    const double* u = _u.data();
    std::vector<double> ret(_grid.n_nodes());
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        ret[i] = _data_by_node[i]->w2(a[i], u[i]);
    }
    return ret;
};
std::vector<double> AssemblerFlux::normalized_area(const std::vector<double>& area) const
{
    std::vector<double> a1(_grid.n_nodes());
    const std::vector<double>& a = (area.size() == 0) ? _area : area;
    for (size_t i = 0; i < _grid.n_nodes(); ++i)
    {
        a1[i] = a[i] / _data[0]->area0 - 1;
    }
    return a1;
}
void AssemblerFlux::reset_flux_calculator(size_t icell, std::shared_ptr<IUpwindFluxCalculator> calc)
{
    _upwind_flux_calculator[icell] = calc;
}

#ifndef ASSEMBLER_FLUX_HPP
#define ASSEMBLER_FLUX_HPP
#include"bflow/fem_grid.hpp"
#include"bflow/flux_calculator.hpp"
#include <memory>
using namespace bflow;

class AssemblerFlux
{
public:
    AssemblerFlux(const FemGrid& grid, const std::vector<const ProblemData*>& data,
                  const std::vector<std::shared_ptr<IUpwindFluxCalculator>>& upwind_calculators);
    AssemblerFlux(const FemGrid& grid, const ProblemData* data,
                  const std::vector<std::shared_ptr<IUpwindFluxCalculator>>& upwind_calculators);
    const CsrMatrix& mass() const;
    const CsrMatrix& block_transport() const;
    void actualize_fluxes(double time, const std::vector<double>& area, const std::vector<double>& velo);
    std::vector<double> dfa_dx() const;
    std::vector<double> dfu_dx() const;
    CsrMatrix viscous_matrix() const;
    CsrMatrix block_u_transport() const;
    std::vector<double> coupling_flux_ua() const;
    std::vector<double> coupling_flux_u2() const;
    std::vector<double> coupling_flux_p() const;
    double compute_residual(const CsrMatrix& lhs, const std::vector<double>& rhs, const std::vector<double>& u) const;
    double vector_norm2(const std::vector<double>& v) const;
    std::vector<double> pressure(const std::vector<double>& area = std::vector<double>()) const;
    std::vector<double> w2() const;
    void reset_flux_calculator(size_t icell, std::shared_ptr<IUpwindFluxCalculator> calc);
private:
    const FemGrid& _grid;
    std::vector<const ProblemData*> _data, _data_by_node;

    double _time = 0;
    CsrMatrix _mass;
    CsrMatrix _tran;
    std::vector<double> _load_vector;
    std::vector<std::shared_ptr<IUpwindFluxCalculator>> _upwind_flux_calculator;
    std::vector<ElementBoundaryFluxes> _upwind_fluxes;
    std::vector<double> _flux_a, _flux_u;
    std::vector<double> _u, _area;
};
#endif

#ifndef __DG_FEM_I_ELEMENT_HPP__
#define __DG_FEM_I_ELEMENT_HPP__

#include "cfd/mat/densemat.hpp"
#include "cfd/grid/i_grid.hpp"
#include "cfd/geom/jacobi.hpp"

namespace dg{

///////////////////////////////////////////////////////////////////////////////
// Element Basis
///////////////////////////////////////////////////////////////////////////////

class IElementBasis{
public:
	virtual ~IElementBasis() = default;

	virtual size_t size() const = 0;
	virtual std::vector<double> parametric_reference_points() const = 0;
	virtual std::vector<double> value(double xi) const = 0;
	virtual std::vector<double> grad(double xi) const = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Element Integrals
///////////////////////////////////////////////////////////////////////////////
class IElementIntegrals{
public:
	virtual ~IElementIntegrals() = default;

	// ⌠
	// ⎮ ϕᵢ dΩ
	// ⌡
	virtual std::vector<double> load_vector() const  { _THROW_NOT_IMP_; }
	// ⌠
	// ⎮ ϕⱼ ϕᵢ dΩ
	// ⌡
	virtual std::vector<double> mass_matrix() const  { _THROW_NOT_IMP_; }
	// ⌠ 𝜕ϕⱼ 𝜕ϕᵢ
	// ⎮ ──  ── dΩ 
	// ⌡ 𝜕x  𝜕x
	virtual std::vector<double> stiff_matrix() const { _THROW_NOT_IMP_; }
	// ⌠ 𝜕ϕⱼ
	// ⎮ ── ϕᵢ dΩ
	// ⌡ 𝜕x
	virtual std::vector<double> dx_matrix() const {_THROW_NOT_IMP_;}
};

///////////////////////////////////////////////////////////////////////////////
// FemElement
///////////////////////////////////////////////////////////////////////////////
struct FemElement{
	double length;
	std::shared_ptr<const IElementBasis> basis;
	std::shared_ptr<const IElementIntegrals> integrals;
};

}
#endif

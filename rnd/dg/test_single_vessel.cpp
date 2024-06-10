#include "common_test.hpp"
#include "mat.hpp"
#include "vtk.hpp"
#include "nonlin_solver.hpp"

namespace{

class FemGrid{
public:
	FemGrid(double length, size_t n, size_t power): _power(power){
		_points.push_back(0);
		_nodes.push_back(0);
		for (size_t i=0; i<n; ++i){
			_points.push_back( double(i+1)/n*length );
			_nodes.push_back(_points.back());
			_nodes.push_back(_points.back());
		}
		_nodes.pop_back();
		if (power > 1){
			_THROW_NOT_IMP_;
		}
	}

	double h() const{
		return _points.back() / n_elements();
	}

	size_t n_nodes() const {
		return _nodes.size();
	}

	size_t n_elements() const {
		return _points.size() - 1;
	}

	size_t n_points() const {
		return _points.size();
	}

	size_t n_local_bases() const {
		return _power + 1;
	}

	double node(size_t i) const{
		return _nodes[i];
	}

	double full_length() const{
		return _points.back() - _points[0];
	}

	size_t closest_node(double x) const {
		double t = 1e16;
		size_t ret = 0;
		for (size_t i=0; i<n_nodes(); ++i){
			double t1 = std::abs(_nodes[i] - x);
			if (t1 < t){
				t = t1;
				ret = i;
			}
		}
		return ret;
	};

	CsrMatrix mass_matrix() const{
		CsrMatrix ret(stencil());
		std::vector<double> local = local_mass_matrix();

		for (size_t ielem=0; ielem<n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_nodes(ielem);

			for (size_t irow=0; irow<n_local_bases(); ++irow)
			for (size_t icol=0; icol<n_local_bases(); ++icol){
				double v = h()/2 * local[irow * n_local_bases() + icol];
				size_t iaddr = ret.get_address(lg[irow], lg[icol]);
				ret.vals()[iaddr] += v;
			}
		}

		return ret;
	}

	CsrMatrix transport_matrix() const{
		CsrMatrix ret(stencil());
		std::vector<double> local = local_transport_matrix();
		for (size_t ielem=0; ielem<n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_nodes(ielem);
			// block diagonal
			for (size_t irow=0; irow<n_local_bases(); ++irow)
			for (size_t icol=0; icol<n_local_bases(); ++icol){
				double v = local[irow * n_local_bases() + icol];
				size_t iaddr = ret.get_address(lg[irow], lg[icol]);
				ret.vals()[iaddr] += v;
			}
		}

		return ret;
	}

	void save_vtk(const std::vector<double>& area,
			const std::vector<double>& velo,
			const std::vector<double>& pressure,
			const std::string s) const {
		save_vtk(s);
		save_vtk_start_point_data(s);
		save_vtk_point_data("area", area, s);
		save_vtk_point_data("velocity", velo, s);
		save_vtk_point_data("pressure", pressure, s);
	}

	void save_vtk(const std::string s) const {
		std::ofstream fs(s);
		fs << "# vtk DataFile Version 3.0" << std::endl;
		fs << "DG" << std::endl;
		fs << "ASCII" << std::endl;
		fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		fs << "POINTS " << n_nodes() << " double" << std::endl;
		for (const double& point: _nodes){
			fs << point << " 0 0" << std::endl;
			fs << std::endl;
		}

		//Cells
		fs << "CELLS  " << n_elements() << "   " << 3 * n_elements() << std::endl;
		for (size_t ielem = 0; ielem < n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_nodes(ielem);
			fs << 2 << " " << lg[0] << " " << lg[1] << std::endl;
		}
		fs << "CELL_TYPES  " << n_elements() << std::endl;
		for (size_t i = 0; i < n_elements(); ++i)
			fs << 3 << std::endl;
		fs.close();
	}
	void save_vtk_start_point_data(const std::string s) const {
		std::ofstream fs(s, std::ios::app);
		fs << "POINT_DATA " << 2*n_elements() << std::endl;
		fs.close();
	}
	void save_vtk_point_data(const std::string data_name, const std::vector<double>& v, const std::string s) const {
		std::ofstream fs(s, std::ios::app);
		fs << "SCALARS " << data_name << " double 1" << std::endl;
		fs << "LOOKUP_TABLE default" << std::endl;
		for (size_t i=0; i<2*n_elements(); ++i){
			fs << v[i] << std::endl;
		}
		fs.close();
	}

	std::vector<size_t> tab_elem_nodes(size_t ielem) const{
		std::vector<size_t> ret {2*ielem, 2*ielem+1};
		if (_power > 1){
			_THROW_NOT_IMP_;
		}
		return ret;
	}

	size_t tab_node_elem(size_t inode) const{
		if (_power > 1){
			_THROW_NOT_IMP_;
		}
		return inode / 2;
	};

	double elem_interpolate(size_t ielem, double x, const std::vector<double>& vec) const {
		std::vector<size_t> n = tab_elem_nodes(ielem);
		double coo = 2 * x / (_nodes[n[1]] - _nodes[n[0]]) - 1;
		if (_power > 1){
			_THROW_NOT_IMP_;
		}
		return (1 - coo)/2 * vec[n[0]] + (1 + coo)/2 * vec[n[1]];
	}
private:
	const size_t _power;
	mutable CsrStencil _stencil;
	std::vector<double> _points;
	std::vector<double> _nodes;

	CsrStencil stencil() const{
		if (_stencil.n_rows() == 0){
			std::vector<std::set<size_t>> lod(n_nodes());
			// block diagonal
			for (size_t ielem=0; ielem < n_elements(); ++ielem){
				std::vector<size_t> lg = tab_elem_nodes(ielem);

				for (size_t irow=0; irow<n_local_bases(); ++irow)
				for (size_t icol=0; icol<n_local_bases(); ++icol){
					lod[lg[irow]].insert(lg[icol]);
				}
			}
			// coupling
			for (size_t ipoint=0; ipoint < n_points() - 1; ++ipoint){
				std::vector<size_t> nodes = tab_point_nodes(ipoint);
				if (nodes.size() == 2){
					lod[nodes[0]].insert(nodes[1]);
					lod[nodes[1]].insert(nodes[0]);
				}
			}
			_stencil.set_stencil(lod);
		}
		return _stencil;
	}

	std::vector<size_t> tab_point_nodes(size_t ipoint) const {
		if (ipoint == 0) return {0};
		else if (ipoint == n_points()-1) return {2*ipoint};
		else return {2*ipoint-1, 2*ipoint};
	}

	std::vector<double> local_mass_matrix() const{
		if (_power == 1){
			return {2.0/3.0, 1.0/3.0, 1.0/3.0, 2.0/3.0};
		} else {
			_THROW_NOT_IMP_;
		}
	}

	std::vector<double> local_transport_matrix() const{
		if (_power == 1){
			return {-0.5, -0.5, 0.5, 0.5};
		} else {
			_THROW_NOT_IMP_;
		}
	}
};

struct ElementBoundaryFluxes{
	double upwind_area_x0 = 0;
	double upwind_area_x1 = 0;
	double upwind_velo_x0 = 0;
	double upwind_velo_x1 = 0;
};

struct ProblemData{
	ProblemData(){
		recompute();
	}

	void recompute(){
		beta = 4.0/3.0*sqrt(pi)*h*E/area0;
		amult = 4*sqrt(beta/2/rho);
		root4_a0 = sqrt(sqrt(area0));
		visc_coef = -2*(profile_order+2)*mu*pi/rho;
	}

	static constexpr double pi = 3.1415926;

	// geometry parameters
	double area0 = pi*1e-4;

	// fluid parameters
	double rho = 1050;
	double mu = 0;
	double profile_order = 9;

	// vessel tissue parameters
	double h = 1.5e-3;
	double E = 4e5;

	// p(a)
	double pressure(double area) const{
		return beta*(std::sqrt(area) - std::sqrt(area0));
	};

	// fluxes
	double flux_a(double area, double vel) const{
		return area*vel;
	};
	double flux_u(double area, double vel) const{
		return vel*vel/2 + pressure(area)/rho;
	}

	// characteristic variables
	double w1(double area, double vel) const{
		return vel + amult*(sqrt(sqrt(area)) - root4_a0);
	}
	double w2(double area, double vel) const{
		return vel - amult*(sqrt(sqrt(area)) - root4_a0);
	}

	// computed characteristics
	double beta;
	double amult;
	double root4_a0;
	double visc_coef;
};

// inflow conditions
double q_inflow1(double t){
	return 1e-6*exp(-1e4*(t-0.05)*(t-0.05));
};
double p_inflow2(double t){
	constexpr double T = 0.33;
	//constexpr double T = 0.01;
	return (T/2 - t > 0) ? 2000*sin(2*ProblemData::pi * t / T) : 0;
};


class IUpwindFluxCalculator{
public:
	virtual void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) = 0;
};

class InflowQFluxCalculator: public IUpwindFluxCalculator{
public:
	InflowQFluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> getq, size_t cell, double eps=1e-12):
		_data(data),
		_getq(getq),
		_cell_right(cell),
		_node_right(grid.tab_elem_nodes(cell)[0]),
		_sys(data.amult, data.root4_a0),
		_eps(eps)
	{
		_a = data.area0;
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double q = _getq();
		double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);
		//double w2_upw = 0;

		double area_upw = _a;
		double velo_upw = q / area_upw;
		_sys.set_qw(q, w2_upw);
		solve_nonlinear_system(_sys, area_upw, velo_upw, _eps);

		fluxes[_cell_right].upwind_area_x0 = area_upw;
		fluxes[_cell_right].upwind_velo_x0 = velo_upw;

		_a = area_upw;
	}
private:
	struct NonlinearSystem: public INonlinearSystem2{
	public:
		NonlinearSystem(double mult, double root4_area0):
			_mult(mult),
			_root4_area0(root4_area0)
		{
			set_qw(0, 0);
		}

		void set_qw(double q, double w2){
			_q = q;
			_w2 = w2;
		};

		std::array<double, 2> f(double area, double velo) const override{
			return {
				area*velo - _q,
				velo - _mult*(sqrt(sqrt(area)) - _root4_area0) - _w2
			};
		};
		std::array<double, 4> jac(double area, double velo) const override{
			return {
				velo, area,
				-0.25*_mult*pow(area, -0.75), 1
			};
		};
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

class InflowPFluxCalculator_NonReflecting: public IUpwindFluxCalculator{
public:
	InflowPFluxCalculator_NonReflecting(const FemGrid& grid, const ProblemData& data, std::function<double()> getp, size_t cell):
		_data(data),
		_getp(getp),
		_cell_right(cell),
		_node_right(grid.tab_elem_nodes(cell)[0])
	{
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double p = _getp();
		//double w2_upw = _getw2();
		double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);
		double a = (p / _data.beta + sqrt(_data.area0));
		a = a*a;
		double w1_upw = 8*sqrt(_data.beta/2/_data.rho)*(sqrt(sqrt(a)) - sqrt(sqrt(_data.area0)));

		double a1 = (w1_upw - w2_upw)/2 / _data.amult + _data.root4_a0;
		double area_upw = a1*a1*a1*a1;
		double velo_upw = (w1_upw + w2_upw)/2.0;

		fluxes[_cell_right].upwind_area_x0 = area_upw;
		fluxes[_cell_right].upwind_velo_x0 = velo_upw;
	}
private:
	const ProblemData& _data;
	std::function<double()> _getp;
	const size_t _cell_right;
	const size_t _node_right;
};

class InflowPFluxCalculator: public IUpwindFluxCalculator{
public:
	InflowPFluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> getp, size_t cell, double eps=1e-12):
		_data(data),
		_getp(getp),
		_cell_right(cell),
		_node_right(grid.tab_elem_nodes(cell)[0]),
		_sys(data.amult, data.area0, data.beta),
		_eps(eps)
	{
		_a = data.area0;
		_v = 0;
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double p = _getp();
		double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);

		double area_upw = _a;
		double velo_upw = _v;
		_sys.set_pw(p, w2_upw);
		solve_nonlinear_system(_sys, area_upw, velo_upw, _eps);

		fluxes[_cell_right].upwind_area_x0 = area_upw;
		fluxes[_cell_right].upwind_velo_x0 = velo_upw;

		_a = area_upw;
		_v = velo_upw;
	}
private:
	struct NonlinearSystem: public INonlinearSystem2{
	public:
		NonlinearSystem(double mult, double area0, double beta):
			_mult(mult),
			_area0(area0),
			_beta(beta),
			_root2_area0(std::pow(area0, 0.5)),
			_root4_area0(std::pow(area0, 0.25))
		{
			set_pw(0, 0);
		}

		void set_pw(double p, double w2){
			_p = p;
			_w2 = w2;
		};

		std::array<double, 2> f(double area, double velo) const override{
			return {
				_beta*(sqrt(area) - _root2_area0) - _p,
				velo - _mult*(sqrt(sqrt(area)) - _root4_area0) - _w2
			};
		};
		std::array<double, 4> jac(double area, double velo) const override{
			return {
				0.5*_beta/sqrt(area), 0,
				-0.25*_mult*pow(area, -0.75), 1
			};
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

class InflowW1FluxCalculator: public IUpwindFluxCalculator{
public:
	InflowW1FluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> getw1, size_t cell):
		_data(data),
		_getw1(getw1),
		_cell_right(cell),
		_node_right(grid.tab_elem_nodes(cell)[0])
	{ }

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double w1_upw = _getw1();
		double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);

		double a1 = (w1_upw - w2_upw)/2 / _data.amult + _data.root4_a0;
		double area_upw = a1*a1*a1*a1;
		double velo_upw = (w1_upw + w2_upw)/2.0;

		fluxes[_cell_right].upwind_area_x0 = area_upw;
		fluxes[_cell_right].upwind_velo_x0 = velo_upw;
	}
private:
	const ProblemData& _data;
	std::function<double()> _getw1;
	const size_t _cell_right;
	const size_t _node_right;
};

class OutflowFluxCalculator: public IUpwindFluxCalculator{
public:
	OutflowFluxCalculator(const FemGrid& grid, const ProblemData& data, size_t cell):
		_data(data),
		_cell_left(cell),
		_node_left(grid.tab_elem_nodes(cell)[1])
	{
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double w1_upw = _data.w1(area[_node_left], velocity[_node_left]);
		double w2_upw = 0;

		double a1 = (w1_upw - w2_upw)/2 / _data.amult + _data.root4_a0;
		double area_upw = a1*a1*a1*a1;
		double velo_upw = (w1_upw + w2_upw)/2.0;

		fluxes[_cell_left].upwind_area_x1 = area_upw;
		fluxes[_cell_left].upwind_velo_x1 = velo_upw;
	}
private:
	const ProblemData& _data;
	const size_t _cell_left;
	const size_t _node_left;
};

class InternalFluxCalculator: public IUpwindFluxCalculator{
public:
	InternalFluxCalculator(const FemGrid& grid, const ProblemData& data, size_t cell_left, size_t cell_right):
		_data(data),
		_node_left(grid.tab_elem_nodes(cell_left)[1]),
		_node_right(grid.tab_elem_nodes(cell_right)[0]),
		_cell_left(cell_left),
		_cell_right(cell_right),
		_amult(data.rho*data.rho/4.0/data.beta/data.beta)
	{
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double w1_upw = _data.w1(area[_node_left], velocity[_node_left]);
		double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);
		double a1 = (w1_upw - w2_upw)/2 / _data.amult + _data.root4_a0;
		double area_upw = a1*a1*a1*a1;
		double velo_upw = (w1_upw + w2_upw)/2.0;

		fluxes[_cell_left].upwind_area_x1 = area_upw;
		fluxes[_cell_left].upwind_velo_x1 = velo_upw;
		fluxes[_cell_right].upwind_area_x0 = area_upw;
		fluxes[_cell_right].upwind_velo_x0 = velo_upw;
	}
private:
	const ProblemData& _data;
	const size_t _node_left;
	const size_t _node_right;
	const size_t _cell_left;
	const size_t _cell_right;
	const double _amult;
};

class MergingFluxCalculator: public IUpwindFluxCalculator{
public:
	MergingFluxCalculator(const FemGrid& grid, const ProblemData& data_left, const ProblemData& data_right, size_t cell_left, size_t cell_right, double eps=1e-12):
		_data_left(data_left),
		_data_right(data_right),
		_node_left(grid.tab_elem_nodes(cell_left)[1]),
		_node_right(grid.tab_elem_nodes(cell_right)[0]),
		_cell_left(cell_left),
		_cell_right(cell_right),
		_sys(data_left, data_right),
		_eps(eps)
	{
		_a1 = _data_left.area0;
		_a2 = _data_right.area0;
		_u1 = 0;
		_u2 = 0;
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double w1_upw = _data_left.w1(area[_node_left], velocity[_node_left]);
		double w2_upw = _data_right.w2(area[_node_right], velocity[_node_right]);

		double area_upw_1 = _a1;
		double velo_upw_1 = _u1;
		double area_upw_2 = _a2;
		double velo_upw_2 = _u2;

		_sys.set_ww(w1_upw, w2_upw);
		solve_nonlinear_system(_sys, area_upw_1, velo_upw_1, area_upw_2, velo_upw_2, _eps);

		//double q = _data_left.flux_a(area_upw_1, velo_upw_1);
		//double p = _data_right.flux_u(area_upw_1, velo_upw_1);
		//velo_upw_2 = q / area_upw_2;

		fluxes[_cell_left].upwind_area_x1 = area_upw_1;
		fluxes[_cell_left].upwind_velo_x1 = velo_upw_1;
		fluxes[_cell_right].upwind_area_x0 = area_upw_2;
		fluxes[_cell_right].upwind_velo_x0 = velo_upw_2;

		_a1 = area_upw_1;
		_u1 = velo_upw_1;
		_a2 = area_upw_2;
		_u2 = velo_upw_2;
	}

private:
	struct NonlinearSystem: public INonlinearSystem4{
	public:
		NonlinearSystem(const ProblemData& data1, const ProblemData& data2):
			_data1(data1), _data2(data2)
		{
			set_ww(0, 0);
		}
		void set_ww(double w1, double w2){
			_w1 = w1;
			_w2 = w2;
		};
		std::array<double, 4> f(double area1, double velo1, double area2, double velo2) const override{
			return {
				_fm*(_data1.flux_a(area1, velo1) - _data2.flux_a(area2, velo2)),
				_fm*(_data1.flux_u(area1, velo1) - _data2.flux_u(area2, velo2)),
				_data1.w1(area1, velo1) - _w1,
				_data2.w2(area2, velo2) - _w2
			};
		};
		std::array<double, 16> jac(double area1, double velo1, double area2, double velo2) const override{
			std::array<double, 16> ret = {
				velo1, area1, -velo2, -area2,
				0.5*_data1.beta/_data1.rho/sqrt(area1), velo1, -0.5*_data2.beta/_data2.rho/sqrt(area2), -velo2,
				sqrt(_data1.beta/2/_data1.rho)*std::pow(area1, -0.75), 1, 0, 0,
				0, 0, -sqrt(_data2.beta/2/_data2.rho)*std::pow(area2, -0.75), 1
			};
			for (size_t i=0; i<8; ++i) ret[i] *= _fm;
			return ret;
		};
	private:
		const ProblemData& _data1;
		const ProblemData& _data2;
		double _w1, _w2;
		static constexpr double _fm = 1;
	};

	const ProblemData& _data_left;
	const ProblemData& _data_right;
	const size_t _node_left, _node_right;
	const size_t _cell_left, _cell_right;
	NonlinearSystem _sys;
	const double _eps;

	double _a1, _a2, _u1, _u2;
};

class Bifurcation3FluxCalculator: public IUpwindFluxCalculator{
public:
	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		_THROW_NOT_IMP_;
	}
};

class Junction3FluxCalculator: public IUpwindFluxCalculator{
public:
	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		_THROW_NOT_IMP_;
	}
};

}

TEST_CASE("Single vessel, inviscid", "[single-vessel-inviscid-explicit]"){
	ProblemData data;
	double time = 0;
	double L = 1;

	FemGrid grid(L, 100, 1);
	double tau = grid.h()/100;
	std::vector<ElementBoundaryFluxes> upwind_fluxes(grid.n_elements());

	VtkUtils::TimeSeriesWriter writer("single-vessel");
	writer.set_time_step(0.005);

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data.area0);
	std::vector<double> pressure(grid.n_nodes(), 0.0);

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
	upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, [&time](){ return q_inflow1(time); }, 0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i-1, i));
	}
	upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements()-1));

	// prepare matrix solver
	CsrMatrix mass = grid.mass_matrix();
	CsrMatrix tran = grid.transport_matrix();
	AmgcMatrixSolver slv;
	slv.set_matrix(mass);

	while (time < 0.2 - 1e-6){
		time += tau;
		//std::cout << "TIME=" << time;
		//std::cout << "  Q=" << data.q_inflow(time) << std::endl;

		// nodewise fluxes
		std::vector<double> flux_a(grid.n_nodes());
		std::vector<double> flux_u(grid.n_nodes());
		for (size_t i=0; i<grid.n_nodes(); ++i){
			flux_a[i] = data.flux_a(area[i], velocity[i]);
			flux_u[i] = data.flux_u(area[i], velocity[i]);
		}

		// upwind fluxes
		for (auto c: upwind_flux_calculator){
			c->compute(area, velocity, upwind_fluxes);
		}

		// assemble rhs
		// E
		std::vector<double> rhs_a = mass.mult_vec(area);
		std::vector<double> rhs_u = mass.mult_vec(velocity);
		// +tau*T
		std::vector<double> tran_a = tran.mult_vec(flux_a);
		std::vector<double> tran_u = tran.mult_vec(flux_u);
		for (size_t i=0; i<grid.n_nodes(); ++i){
			rhs_a[i] += tau * tran_a[i];
			rhs_u[i] += tau * tran_u[i];
		}
		// + coupling
		for (size_t ielem=0; ielem<grid.n_elements(); ++ielem){
			size_t node0 = grid.tab_elem_nodes(ielem)[0];
			size_t node1 = grid.tab_elem_nodes(ielem)[1];
			const auto& uf = upwind_fluxes[ielem];
			rhs_a[node0] += tau*uf.upwind_area_x0 * uf.upwind_velo_x0;
			rhs_a[node1] -= tau*uf.upwind_area_x1 * uf.upwind_velo_x1;
			rhs_u[node0] += tau*(0.5 * uf.upwind_velo_x0 * uf.upwind_velo_x0 + data.pressure(uf.upwind_area_x0)/data.rho);
			rhs_u[node1] -= tau*(0.5 * uf.upwind_velo_x1 * uf.upwind_velo_x1 + data.pressure(uf.upwind_area_x1)/data.rho);
		}
		//solve
		slv.solve(rhs_a, area);
		slv.solve(rhs_u, velocity);

		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) {
			grid.save_vtk(area, velocity, pressure, out_filename);
		}
	}

	double maxp = *std::max_element(pressure.begin(), pressure.end());
	CHECK(maxp == Approx(24.3958).margin(1e-3));
}

namespace{

class Assembler{
public:
	Assembler(
		const FemGrid& grid,
		const std::vector<const ProblemData*>& data,
		const std::vector<std::shared_ptr<IUpwindFluxCalculator>>& upwind_calculators)
		: _grid(grid), _data(data), _upwind_flux_calculator(upwind_calculators)
	{
		_data_by_node = std::vector<const ProblemData*>(grid.n_nodes());
		for (size_t i=0; i<grid.n_nodes(); ++i){
			size_t icell = grid.tab_node_elem(i);
			_data_by_node[i] = _data[icell];
		}

		_mass = _grid.mass_matrix();
		_tran = _grid.transport_matrix();

		_upwind_fluxes.resize(grid.n_elements());
		_flux_a.resize(_grid.n_nodes());
		_flux_u.resize(_grid.n_nodes());

		_load_vector = _mass.mult_vec(std::vector<double>(_grid.n_nodes(), 1.0));
	}

	Assembler(
		const FemGrid& grid,
		const ProblemData* data,
		const std::vector<std::shared_ptr<IUpwindFluxCalculator>>& upwind_calculators)
		: Assembler(grid, std::vector<const ProblemData*>(grid.n_elements(), data), upwind_calculators){}

	const CsrMatrix& mass() const{
		return _mass;
	}

	const CsrMatrix& block_transport() const{
		return _tran;
	}

	void actualize_fluxes(double time, const std::vector<double>& area, const std::vector<double>& velo){
		_time = time;
		_u = velo;
		_area = area;

		// compute upwind fluxes
		for (auto c: _upwind_flux_calculator){
			c->compute(area, velo, _upwind_fluxes);
		}

		// nodewise fluxes
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			_flux_a[i] = _data_by_node[i]->flux_a(area[i], velo[i]);
			_flux_u[i] = _data_by_node[i]->flux_u(area[i], velo[i]);
		}
	}

	std::vector<double> dfa_dx() const{
		// nodewise flux
		std::vector<double> ret = _tran.mult_vec(_flux_a);
		for (auto& v: ret) v*=-1;

		// coupling
		for (size_t ielem=0; ielem<_grid.n_elements(); ++ielem){
			size_t node0 = _grid.tab_elem_nodes(ielem)[0];
			size_t node1 = _grid.tab_elem_nodes(ielem)[1];
			const auto& uf = _upwind_fluxes[ielem];
			ret[node0] -= uf.upwind_area_x0 * uf.upwind_velo_x0;
			ret[node1] += uf.upwind_area_x1 * uf.upwind_velo_x1;
		}

		return ret;
	}
	std::vector<double> dfu_dx() const{
		// nodewise flux
		std::vector<double> ret = _tran.mult_vec(_flux_u);
		for (auto& v: ret) v*=-1;

		// coupling
		for (size_t ielem=0; ielem<_grid.n_elements(); ++ielem){
			size_t node0 = _grid.tab_elem_nodes(ielem)[0];
			size_t node1 = _grid.tab_elem_nodes(ielem)[1];
			const auto& uf = _upwind_fluxes[ielem];

			ret[node0] -= 0.5 * uf.upwind_velo_x0 * uf.upwind_velo_x0 + _data[ielem]->pressure(uf.upwind_area_x0)/_data[ielem]->rho;
			ret[node1] += 0.5 * uf.upwind_velo_x1 * uf.upwind_velo_x1 + _data[ielem]->pressure(uf.upwind_area_x1)/_data[ielem]->rho;
		}
		return ret;
	}

	CsrMatrix viscous_matrix() const{
		CsrMatrix ret(_mass);
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			for (size_t a=ret.addr()[i]; a < ret.addr()[i+1]; ++a){
				size_t col = ret.cols()[a];
				ret.vals()[a] *= (_data_by_node[i]->visc_coef / _area[col]);
			}
		}
		return ret;
	}

	CsrMatrix block_u_transport() const {
		CsrMatrix ret(_tran);
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			for (size_t a=ret.addr()[i]; a < ret.addr()[i+1]; ++a){
				size_t col = ret.cols()[a];
				ret.vals()[a] *= _u[col];
			}
		}
		return ret;
	}

	std::vector<double> coupling_flux_ua() const{
		std::vector<double> ret(_grid.n_nodes(), 0.0);
		for (size_t ielem=0; ielem<_grid.n_elements(); ++ielem){
			size_t node0 = _grid.tab_elem_nodes(ielem)[0];
			size_t node1 = _grid.tab_elem_nodes(ielem)[1];
			ret[node0] -= _upwind_fluxes[ielem].upwind_area_x0 * _upwind_fluxes[ielem].upwind_velo_x0;
			ret[node1] += _upwind_fluxes[ielem].upwind_area_x1 * _upwind_fluxes[ielem].upwind_velo_x1;
		}
		return ret;
	}

	std::vector<double> coupling_flux_u2() const{
		std::vector<double> ret(_grid.n_nodes(), 0.0);
		for (size_t ielem=0; ielem<_grid.n_elements(); ++ielem){
			size_t node0 = _grid.tab_elem_nodes(ielem)[0];
			size_t node1 = _grid.tab_elem_nodes(ielem)[1];
			ret[node0] -= _upwind_fluxes[ielem].upwind_velo_x0 * _upwind_fluxes[ielem].upwind_velo_x0;
			ret[node1] += _upwind_fluxes[ielem].upwind_velo_x1 * _upwind_fluxes[ielem].upwind_velo_x1;
		}
		return ret;
	}

	std::vector<double> coupling_flux_p() const{
		std::vector<double> ret(_grid.n_nodes(), 0.0);
		for (size_t ielem=0; ielem<_grid.n_elements(); ++ielem){
			size_t node0 = _grid.tab_elem_nodes(ielem)[0];
			size_t node1 = _grid.tab_elem_nodes(ielem)[1];
			ret[node0] -= _data[ielem]->pressure(_upwind_fluxes[ielem].upwind_area_x0);
			ret[node1] += _data[ielem]->pressure(_upwind_fluxes[ielem].upwind_area_x1);
		}
		return ret;
	}

	double compute_residual(const CsrMatrix& lhs, const std::vector<double>& rhs, const std::vector<double>& u) const{
		std::vector<double> r = lhs.mult_vec(u);
		for (size_t i=0; i<r.size(); ++i){
			r[i] = rhs[i] - r[i];
		}
		return vector_norm2(r);
	}

	double vector_norm2(const std::vector<double>& v) const{
		double s = 0;
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			s += v[i] * v[i] * _load_vector[i];
		}
		return std::sqrt(s/_grid.full_length());
	}

	std::vector<double> pressure(const std::vector<double>& area=std::vector<double>()) const{
		const double* a = (area.size() == 0) ? _area.data() : area.data();
		std::vector<double> ret(_grid.n_nodes());
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			ret[i] = _data_by_node[i]->pressure(a[i]);
		}
		return ret;
	};

	std::vector<double> w2() const{
		const double* a = _area.data();
		const double* u = _u.data();
		std::vector<double> ret(_grid.n_nodes());
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			ret[i] = _data_by_node[i]->w2(a[i], u[i]);
		}
		return ret;
	};

	void save_vtk(std::string filename,
			const std::vector<double>& area=std::vector<double>(),
			const std::vector<double>& velo=std::vector<double>()){

		const std::vector<double>& a = (area.size() == 0) ? _area : area;
		const std::vector<double>& u = (velo.size() == 0) ? _u : velo;

		std::vector<double> p = pressure(a);
		std::vector<double> w1(_grid.n_nodes()), w2(_grid.n_nodes()), a1(_grid.n_nodes());
		for (size_t i=0; i<_grid.n_nodes(); ++i){
			w1[i] = _data_by_node[i]->w1(a[i], u[i]);
			w2[i] = _data_by_node[i]->w2(a[i], u[i]);
			a1[i] = a[i] / _data[0]->area0 - 1;
		}

		_grid.save_vtk(filename);
		_grid.save_vtk_start_point_data(filename);
		_grid.save_vtk_point_data("area", a1, filename);
		_grid.save_vtk_point_data("pressure", p, filename);
		_grid.save_vtk_point_data("velocity", u, filename);
		_grid.save_vtk_point_data("w1", w1, filename);
		_grid.save_vtk_point_data("w2", w2, filename);
	};

	void reset_flux_calculator(size_t icell, std::shared_ptr<IUpwindFluxCalculator> calc){
		_upwind_flux_calculator[icell] = calc;
	}
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

};

TEST_CASE("Single vessel, inviscid, ver2", "[single-vessel-inviscid-explicit2]"){
	ProblemData data;
	double L = 1;
	data.recompute();
	double time = 0;

	FemGrid grid(L, L*100, 1);
	double tau = grid.h()/100;

	VtkUtils::TimeSeriesWriter writer("single-vessel");
	writer.set_time_step(0.005);

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data.area0);
	std::vector<double> pressure(grid.n_nodes(), 0.0);

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
	upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, [&time](){ return q_inflow1(time); }, 0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i-1, i));
	}
	upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements()-1));

	Assembler assem(grid, &data, upwind_flux_calculator);
	CHECK(assem.vector_norm2(std::vector<double>(grid.n_nodes(), 3)) == Approx(3.0).margin(1e-12));

	// prepare matrix solver
	AmgcMatrixSolver slv;
	slv.set_matrix(assem.mass());

	while (time < 0.2 - 1e-6){
		assem.actualize_fluxes(time, area, velocity);
		std::vector<double> dfadx = assem.dfa_dx();
		std::vector<double> dfudx = assem.dfu_dx();

		time += tau;
		//std::cout << "TIME=" << time;
		//std::cout << "  Q=" << data.q_inflow(time) << std::endl;

		// rhs_a  = A - tau*d(uA)/dx
		std::vector<double> rhs_a = assem.mass().mult_vec(area);
		for (size_t i=0; i<rhs_a.size(); ++i){
			rhs_a[i] -= tau*dfadx[i];
		}

		// rhs1_u = u
		std::vector<double> rhs_u = assem.mass().mult_vec(velocity);
		//          -tau*(1-theta)d(Pt)/dx
		for (size_t i=0; i<rhs_u.size(); ++i){
			rhs_u[i] -= tau*dfudx[i];
		}
		//          +tau*(1-theta)*kr*u/a
		//for (size_t i=0; i<rhs_u.size(); ++i){
		//        rhs_u[i] += tau*visc_term[i];
		//}

		slv.solve(rhs_a, area);
		slv.solve(rhs_u, velocity);

		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) {
			grid.save_vtk(area, velocity, pressure, out_filename);
		}
	}

	double maxp = *std::max_element(pressure.begin(), pressure.end());
	CHECK(maxp == Approx(24.3632).margin(1e-3));
}


TEST_CASE("Single vessel, inviscid, implicit", "[single-vessel-inviscid-implicit]"){
	ProblemData data;
	double L = 1;
	data.recompute();
	double time = 0;

	FemGrid grid(L, L*30, 1);
	double tau = grid.h()/100;


	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data.area0);
	std::vector<double> pressure(grid.n_nodes(), 0.0);
	for (size_t i=0; i<grid.n_nodes(); ++i){
		pressure[i] = data.pressure(area[i]);
	}

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
	upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, [&time](){ return q_inflow1(time); }, 0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i-1, i));
	}
	upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements()-1));

	Assembler assem(grid, &data, upwind_flux_calculator);

	AmgcMatrixSolver slv;

	VtkUtils::TimeSeriesWriter writer("single-vessel");
	writer.set_time_step(0.005);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) {
		grid.save_vtk(area, velocity, pressure, out_filename);
	}

	size_t iter_max = 100;

	while (time < 0.2 - 1e-6){
		// assemble right hand side
		std::vector<double> rhs1_a = assem.mass().mult_vec(area);
		std::vector<double> rhs1_u = assem.mass().mult_vec(velocity);

		time += tau;
		//std::cout << "TIME=" << time;
		//std::cout << "  Q=" << data.q_inflow(time) << std::endl;

		double err_a = 1e6;
		double err_u = 1e6;
		double err_max = 1e-10;

		for (size_t it=0; it < iter_max; ++it){
			// 1. --- Area equation
			assem.actualize_fluxes(time, area, velocity);
			// 1.1 LHS
			CsrMatrix lhs_a = assem.mass();
			CsrMatrix block_u_transport = assem.block_u_transport();
			for (size_t i=0; i<lhs_a.n_nonzeros(); ++i){
				lhs_a.vals()[i] -= tau * block_u_transport.vals()[i];
			}

			// 1.2 RHS
			std::vector<double> rhs_a(rhs1_a);
			// coupling flux
			std::vector<double> coupling_flux_ua = assem.coupling_flux_ua();
			for (size_t i=0; i<rhs_a.size(); ++i){
				rhs_a[i] -= tau*coupling_flux_ua[i];
			}
			// 1.3 residual
			if (it > 0){
				err_a = assem.compute_residual(lhs_a, rhs_a, area);
			}

			// 1.4 Solve
			slv.set_matrix(lhs_a);
			slv.solve(rhs_a, area);

			// 2. ---- Velocity equation
			assem.actualize_fluxes(time, area, velocity);
			// 2.1 LHS
			CsrMatrix lhs_u = assem.mass();
			block_u_transport = assem.block_u_transport();
			for (size_t i=0; i<lhs_u.n_nonzeros(); ++i){
				lhs_u.vals()[i] -= 0.5 * tau * block_u_transport.vals()[i];
			}
			// 2.2 RHS
			std::vector<double> rhs_u(rhs1_u);
			// coupling flux
			std::vector<double> coupling_flux_u2 = assem.coupling_flux_u2();
			std::vector<double> coupling_flux_p = assem.coupling_flux_p();
			for (size_t i=0; i<rhs_a.size(); ++i){
				rhs_u[i] -= tau*(0.5*coupling_flux_u2[i] + coupling_flux_p[i]/data.rho);
			}
			// block p
			std::vector<double> p(grid.n_nodes());
			for (size_t i=0; i<grid.n_nodes(); ++i){
				p[i] = data.pressure(area[i]);
			}
			p = assem.block_transport().mult_vec(p);
			for (size_t i=0; i<rhs_a.size(); ++i){
				rhs_u[i] += tau * p[i] / data.rho;
			}

			// 2.3 residual
			if (it > 0){
				err_u = assem.compute_residual(lhs_u, rhs_u, velocity);
			}
			// 2.4 Solve
			slv.set_matrix(lhs_u);
			slv.solve(rhs_u, velocity);

			// break conditions
			if (it > 0){
				//std::cout << err_u << " " << err_a << std::endl;
				if (std::max(err_u, err_a) < err_max){
					//std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				} else if (it == iter_max - 1){
					std::cout << "Warning: internal iterations did not converge: " << err_u << ", " << err_a << std::endl;
				}
			}
		}

		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) {
			grid.save_vtk(area, velocity, pressure, out_filename);
		}
	}

	double maxp = *std::max_element(pressure.begin(), pressure.end());
	std::cout << maxp << std::endl;
	CHECK(maxp == Approx(14.0287775452).margin(1e-3));
}


TEST_CASE("Single vessel, inviscid, theta-scheme", "[single-vessel-inviscid-theta]"){
	ProblemData data;
	double L = 1;
	data.recompute();
	double time = 0;
	double theta = 0.5;

	FemGrid grid(L, L*30, 1);
	double tau = grid.h()/100;

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data.area0);
	std::vector<double> pressure(grid.n_nodes(), 0.0);
	for (size_t i=0; i<grid.n_nodes(); ++i){
		pressure[i] = data.pressure(area[i]);
	}

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
	upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, [&time](){ return q_inflow1(time); }, 0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i-1, i));
	}
	upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements()-1));

	Assembler assem(grid, &data, upwind_flux_calculator);

	AmgcMatrixSolver slv;

	VtkUtils::TimeSeriesWriter writer("single-vessel");
	writer.set_time_step(0.005);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) {
		grid.save_vtk(area, velocity, pressure, out_filename);
	}

	size_t iter_max = 10;

	CsrMatrix block_u_transport;
	std::vector<double> coupling_flux_ua;
	std::vector<double> coupling_flux_u2;
	std::vector<double> coupling_flux_p;
	std::vector<double> tmp;

	while (time < 0.2 - 1e-6){
		// assemble right hand side
		assem.actualize_fluxes(time, area, velocity);
		block_u_transport = assem.block_u_transport();
		// 1. ---- Area equation
		std::vector<double> rhs1_a = assem.mass().mult_vec(area);
		for (size_t i=0; i<grid.n_nodes(); ++i){
			rhs1_a[i] += (1 - theta) * tau * block_u_transport.mult_vec(i, area);
		}
		coupling_flux_ua = assem.coupling_flux_ua();
		for (size_t i=0; i<rhs1_a.size(); ++i){
			rhs1_a[i] -= (1 - theta) * tau * coupling_flux_ua[i];
		}
		// 2. ---- Velocity equation
		std::vector<double> rhs1_u = assem.mass().mult_vec(velocity);
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] += 0.5 * (1 - theta) * tau * block_u_transport.mult_vec(i, velocity);
		}
		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		tmp = assem.block_transport().mult_vec(pressure);
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] += (1 - theta) * tau * tmp[i] / data.rho;
		}
		coupling_flux_u2 = assem.coupling_flux_u2();
		coupling_flux_p = assem.coupling_flux_p();
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] -= (1 - theta) * tau * (coupling_flux_u2[i]/2 + coupling_flux_p[i]/data.rho);
		}

		time += tau;
		//std::cout << "TIME=" << time;
		//std::cout << "  Q=" << data.q_inflow(time) << std::endl;

		double err_a = 1e6;
		double err_u = 1e6;
		double err_max = 1e-10;

		for (size_t it=0; it < iter_max; ++it){
			// 1. --- Area equation
			assem.actualize_fluxes(time, area, velocity);
			// 1.1 LHS
			CsrMatrix lhs_a = assem.mass();
			lhs_a.plus(-theta * tau, assem.block_u_transport());
			// 1.2 RHS
			std::vector<double> rhs_a(rhs1_a);
			// coupling flux
			coupling_flux_ua = assem.coupling_flux_ua();
			for (size_t i=0; i<rhs_a.size(); ++i){
				rhs_a[i] -= theta * tau*coupling_flux_ua[i];
			}
			// 1.3 residual
			if (it > 0){
				err_a = assem.compute_residual(lhs_a, rhs_a, area);
			}

			// 1.4 Solve
			slv.set_matrix(lhs_a);
			slv.solve(rhs_a, area);

			// 2. ---- Velocity equation
			assem.actualize_fluxes(time, area, velocity);
			// 2.1 LHS
			CsrMatrix lhs_u = assem.mass();
			lhs_u.plus(-0.5 * theta * tau, assem.block_u_transport());
			// 2.2 RHS
			std::vector<double> rhs_u(rhs1_u);
			// coupling flux
			coupling_flux_u2 = assem.coupling_flux_u2();
			coupling_flux_p = assem.coupling_flux_p();
			for (size_t i=0; i<rhs_u.size(); ++i){
				rhs_u[i] -= theta * tau*(0.5*coupling_flux_u2[i] + coupling_flux_p[i]/data.rho);
			}
			// block p
			for (size_t i=0; i<grid.n_nodes(); ++i){
				pressure[i] = data.pressure(area[i]);
			}
			tmp = assem.block_transport().mult_vec(pressure);
			for (size_t i=0; i<rhs_u.size(); ++i){
				rhs_u[i] += theta * tau * tmp[i] / data.rho;
			}

			// 2.3 residual
			if (it > 0){
				err_u = assem.compute_residual(lhs_u, rhs_u, velocity);
			}
			// 2.4 Solve
			slv.set_matrix(lhs_u);
			slv.solve(rhs_u, velocity);

			// break conditions
			if (it > 0){
				//std::cout << err_u << " " << err_a << std::endl;
				if (std::max(err_u, err_a) < err_max){
					//std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				} else if (it == iter_max - 1){
					std::cout << "Warning: internal iterations did not converge: " << err_u << ", " << err_a << std::endl;
				}
			}
		}
		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) {
			grid.save_vtk(area, velocity, pressure, out_filename);
		}
	}
	double maxp = *std::max_element(pressure.begin(), pressure.end());
	CHECK(maxp == Approx(18.117256587).margin(1e-3));
}

namespace {
void time_value_vtk(double tau, const std::vector<double>& v, std::string filename){
	std::ofstream fs(filename);
	size_t n_nodes = v.size();
	size_t n_cells = v.size()-1;
	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << "DG" << std::endl;
	fs << "ASCII" << std::endl;
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << n_nodes << " double" << std::endl;
	for (size_t i=0; i<n_nodes; ++i){
		fs << i*tau << " 0 0" << std::endl;
	}

	//Cells
	fs << "CELLS  " << n_cells << "   " << 3 * n_cells << std::endl;
	for (size_t ielem = 0; ielem < n_cells; ++ielem){
		fs << 2 << " " << ielem << " " << ielem + 1 << std::endl;
	}
	fs << "CELL_TYPES  " << n_cells << std::endl;
	for (size_t i = 0; i < n_cells; ++i)
		fs << 3 << std::endl;

	// Data
	fs << "POINT_DATA " << n_nodes << std::endl;
	fs << "SCALARS area  double 1" << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	for (size_t i=0; i<n_nodes; ++i){
		fs << v[i] << std::endl;
	}
	fs.close();
}
}

TEST_CASE("Single vessel, different beta properties", "[2props-vessel]"){
	ProblemData data1;
	data1.area0 = ProblemData::pi / 4.0;
	data1.rho = 1;
	data1.h = 1;
	data1.E = 84'628.5 * std::sqrt(ProblemData::pi);
	data1.recompute();

	ProblemData data2;
	data2.area0 = ProblemData::pi / 4.0;
	data2.rho = 1;
	data2.h = 1;
	data2.E = data1.E * 100;
	data2.recompute();

	double time = 0;
	double theta = 0.5;
	double L = 15;
	size_t k3 = 10;

	FemGrid grid(L, 3*k3, 1);
	std::vector<int> cell_types(grid.n_elements(), 1);
	for (size_t i=k3; i<2*k3; ++i){
		cell_types[i] = 2;
	}
	double tau = grid.h()/50000;

	size_t monitoring_node1 = grid.closest_node(0.25 * L);
	size_t monitoring_node2 = grid.closest_node(0.5 * L);
	size_t monitoring_node3 = grid.closest_node(0.75 * L);
	std::vector<double> monitor1, monitor2, monitor3;

	std::vector<const ProblemData*> data(grid.n_elements());
	for (size_t icell=0; icell<grid.n_elements(); ++icell){
		data[icell] =  (cell_types[icell] == 1) ? &data1 : &data2;
	}

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> flux_calculators(grid.n_points());
	flux_calculators[0].reset(
		new InflowPFluxCalculator_NonReflecting(grid, data1,
			[&time](){ return p_inflow2(time); },
			0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		size_t cell_left = i-1;
		size_t cell_right = i;
		ProblemData* data_left = (cell_types[cell_left] == 1) ? &data1 : &data2;
		ProblemData* data_right = (cell_types[cell_right] == 1) ? &data1 : &data2;
		if (data_left == data_right){
			flux_calculators[i].reset(new InternalFluxCalculator(grid, *data_left, cell_left, cell_right));
		} else {
			flux_calculators[i].reset(new MergingFluxCalculator(grid, *data_left, *data_right, cell_left, cell_right, 1e-8));
		}
	}
	flux_calculators.back().reset(new OutflowFluxCalculator(grid, data1, grid.n_elements()-1));
	Assembler assem(grid, data, flux_calculators);

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data1.area0);

	AmgcMatrixSolver slv;

	VtkUtils::TimeSeriesWriter writer("single-vessel");
	writer.set_time_step(0.005);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) {
		assem.save_vtk(out_filename, area, velocity);
	}

	size_t iter_max = 1;

	CsrMatrix block_u_transport;
	std::vector<double> coupling_flux_ua;
	std::vector<double> coupling_flux_u2;
	std::vector<double> coupling_flux_p;
	std::vector<double> tmp;


	//while (time < 0.25 - 1e-12){
	while (time < 0.05 - 1e-12){
		// assemble right hand side
		assem.actualize_fluxes(time, area, velocity);
		block_u_transport = assem.block_u_transport();
		// 1. ---- Area equation
		std::vector<double> rhs1_a = assem.mass().mult_vec(area);
		for (size_t i=0; i<grid.n_nodes(); ++i){
			rhs1_a[i] += (1 - theta) * tau * block_u_transport.mult_vec(i, area);
		}
		coupling_flux_ua = assem.coupling_flux_ua();
		for (size_t i=0; i<rhs1_a.size(); ++i){
			rhs1_a[i] -= (1 - theta) * tau * coupling_flux_ua[i];
		}
		// 2. ---- Velocity equation
		std::vector<double> rhs1_u = assem.mass().mult_vec(velocity);
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] += 0.5 * (1 - theta) * tau * block_u_transport.mult_vec(i, velocity);
		}
		tmp = assem.block_transport().mult_vec(assem.pressure());
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] += (1 - theta) * tau * tmp[i] / data1.rho;
		}
		coupling_flux_u2 = assem.coupling_flux_u2();
		coupling_flux_p = assem.coupling_flux_p();
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] -= (1 - theta) * tau * (coupling_flux_u2[i]/2 + coupling_flux_p[i]/data1.rho);
		}

		time += tau;
		//std::cout << "TIME=" << time;
		//std::cout << "  P=" << p_inflow2(time) << std::endl;

		double err_a = 1e6;
		double err_u = 1e6;
		double err_max = 1e-8;

		for (size_t it=0; it < iter_max; ++it){
			// 1. --- Area equation
			assem.actualize_fluxes(time, area, velocity);
			// 1.1 LHS
			CsrMatrix lhs_a = assem.mass();
			lhs_a.plus(-theta * tau, assem.block_u_transport());
			// 1.2 RHS
			std::vector<double> rhs_a(rhs1_a);
			// coupling flux
			coupling_flux_ua = assem.coupling_flux_ua();
			for (size_t i=0; i<rhs_a.size(); ++i){
				rhs_a[i] -= theta * tau*coupling_flux_ua[i];
			}

			// 1.3 residual
			if (it > 0){
				err_a = assem.compute_residual(lhs_a, rhs_a, area);
			}

			// 1.4 Solve
			slv.set_matrix(lhs_a);
			slv.solve(rhs_a, area);

			// 2. ---- Velocity equation
			assem.actualize_fluxes(time, area, velocity);
			// 2.1 LHS
			CsrMatrix lhs_u = assem.mass();
			lhs_u.plus(-0.5 * theta * tau, assem.block_u_transport());
			// 2.2 RHS
			std::vector<double> rhs_u(rhs1_u);
			// coupling flux
			coupling_flux_u2 = assem.coupling_flux_u2();
			coupling_flux_p = assem.coupling_flux_p();
			for (size_t i=0; i<rhs_u.size(); ++i){
				rhs_u[i] -= theta * tau*(0.5*coupling_flux_u2[i] + coupling_flux_p[i]/data1.rho);
			}
			// block p
			tmp = assem.block_transport().mult_vec(assem.pressure());
			for (size_t i=0; i<rhs_u.size(); ++i){
				rhs_u[i] += theta * tau * tmp[i] / data1.rho;
			}

			// 2.3 residual
			if (it > 0){
				err_u = assem.compute_residual(lhs_u, rhs_u, velocity);
			}
			// 2.4 Solve
			slv.set_matrix(lhs_u);
			slv.solve(rhs_u, velocity);

			// break conditions
			if (it > 0){
				//std::cout << err_u << " " << err_a << std::endl;
				if (std::max(err_u, err_a) < err_max){
					//std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				} else if (it == iter_max - 1){
					std::cout << "Warning: internal iterations did not converge: " << err_u << ", " << err_a << std::endl;
				}
			}
		}
		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) {
			assem.save_vtk(out_filename, area, velocity);
		}

		std::vector<double> pressure = assem.pressure(area);
		monitor1.push_back(pressure[monitoring_node1]);
		monitor2.push_back(pressure[monitoring_node2]);
		monitor3.push_back(pressure[monitoring_node3]);
	}

	std::vector<double> pressure = assem.pressure(area);
	CHECK(pressure[20] == Approx(1521.4710994706).margin(1e-3));

	time_value_vtk(tau, monitor1, "monitor1.vtk");
	time_value_vtk(tau, monitor2, "monitor2.vtk");
	time_value_vtk(tau, monitor3, "monitor3.vtk");
}


TEST_CASE("Single vessel, viscous", "[single-vessel-viscous]"){
	ProblemData data;
	data.mu = 0.004;
	data.recompute();

	double L = 5;
	double time = 0;
	double theta = 0.5;

	FemGrid grid(L, L*30, 1);
	double tau = grid.h()/50;

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data.area0);
	std::vector<double> pressure(grid.n_nodes(), 0.0);
	for (size_t i=0; i<grid.n_nodes(); ++i){
		pressure[i] = data.pressure(area[i]);
	}

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
	upwind_flux_calculator[0].reset(new InflowQFluxCalculator(grid, data, [&time](){ return q_inflow1(time); }, 0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i-1, i));
	}
	upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements()-1));

	Assembler assem(grid, &data, upwind_flux_calculator);

	AmgcMatrixSolver slv;

	VtkUtils::TimeSeriesWriter writer("single-vessel");
	writer.set_time_step(0.05);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) {
		grid.save_vtk(area, velocity, pressure, out_filename);
	}

	size_t iter_max = 10;

	CsrMatrix block_u_transport;
	std::vector<double> coupling_flux_ua;
	std::vector<double> coupling_flux_u2;
	std::vector<double> coupling_flux_p;
	std::vector<double> tmp;

	while (time < 0.3 - 1e-6){
		// assemble right hand side
		assem.actualize_fluxes(time, area, velocity);
		block_u_transport = assem.block_u_transport();
		// 1. ---- Area equation
		std::vector<double> rhs1_a = assem.mass().mult_vec(area);
		for (size_t i=0; i<grid.n_nodes(); ++i){
			rhs1_a[i] += (1 - theta) * tau * block_u_transport.mult_vec(i, area);
		}
		coupling_flux_ua = assem.coupling_flux_ua();
		for (size_t i=0; i<rhs1_a.size(); ++i){
			rhs1_a[i] -= (1 - theta) * tau * coupling_flux_ua[i];
		}
		// 2. ---- Velocity equation
		std::vector<double> rhs1_u = assem.mass().mult_vec(velocity);
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] += 0.5 * (1 - theta) * tau * block_u_transport.mult_vec(i, velocity);
		}
		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		tmp = assem.block_transport().mult_vec(pressure);
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] += (1 - theta) * tau * tmp[i] / data.rho;
		}
		coupling_flux_u2 = assem.coupling_flux_u2();
		coupling_flux_p = assem.coupling_flux_p();
		for (size_t i=0; i<rhs1_u.size(); ++i){
			rhs1_u[i] -= (1 - theta) * tau * (coupling_flux_u2[i]/2 + coupling_flux_p[i]/data.rho);
		}
		tmp = assem.viscous_matrix().mult_vec(velocity);
		for (size_t i=0; i<grid.n_nodes(); ++i){
			rhs1_u[i] += (1-theta) * tau * tmp[i];
		}

		time += tau;
		//std::cout << "TIME=" << time;
		//std::cout << "  Q=" << q_inflow1(time) << std::endl;

		double err_a = 1e6;
		double err_u = 1e6;
		double err_max = 1e-10;

		for (size_t it=0; it < iter_max; ++it){
			// 1. --- Area equation
			assem.actualize_fluxes(time, area, velocity);
			// 1.1 LHS
			CsrMatrix lhs_a = assem.mass();
			lhs_a.plus(-theta * tau, assem.block_u_transport());
			// 1.2 RHS
			std::vector<double> rhs_a(rhs1_a);
			// coupling flux
			coupling_flux_ua = assem.coupling_flux_ua();
			for (size_t i=0; i<rhs_a.size(); ++i){
				rhs_a[i] -= theta * tau*coupling_flux_ua[i];
			}
			// 1.3 residual
			if (it > 0){
				err_a = assem.compute_residual(lhs_a, rhs_a, area);
			}

			// 1.4 Solve
			slv.set_matrix(lhs_a);
			slv.solve(rhs_a, area);

			// 2. ---- Velocity equation
			assem.actualize_fluxes(time, area, velocity);
			// 2.1 LHS
			CsrMatrix lhs_u = assem.mass();
			lhs_u.plus(-0.5 * theta * tau, assem.block_u_transport());
			lhs_u.plus(-theta * tau, assem.viscous_matrix());
			// 2.2 RHS
			std::vector<double> rhs_u(rhs1_u);
			// coupling flux
			coupling_flux_u2 = assem.coupling_flux_u2();
			coupling_flux_p = assem.coupling_flux_p();
			for (size_t i=0; i<rhs_u.size(); ++i){
				rhs_u[i] -= theta * tau*(0.5*coupling_flux_u2[i] + coupling_flux_p[i]/data.rho);
			}
			// block p
			for (size_t i=0; i<grid.n_nodes(); ++i){
				pressure[i] = data.pressure(area[i]);
			}
			tmp = assem.block_transport().mult_vec(pressure);
			for (size_t i=0; i<rhs_u.size(); ++i){
				rhs_u[i] += theta * tau * tmp[i] / data.rho;
			}

			// 2.3 residual
			if (it > 0){
				err_u = assem.compute_residual(lhs_u, rhs_u, velocity);
			}
			// 2.4 Solve
			slv.set_matrix(lhs_u);
			slv.solve(rhs_u, velocity);

			// break conditions
			if (it > 0){
				//std::cout << err_u << " " << err_a << std::endl;
				if (std::max(err_u, err_a) < err_max){
					//std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				} else if (it == iter_max - 1){
					std::cout << "Warning: internal iterations did not converge: " << err_u << ", " << err_a << std::endl;
				}
			}
		}
		for (size_t i=0; i<grid.n_nodes(); ++i){
			pressure[i] = data.pressure(area[i]);
		}
		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) {
			grid.save_vtk(area, velocity, pressure, out_filename);
		}
	}
	double maxp = *std::max_element(pressure.begin(), pressure.end());
	CHECK(maxp == Approx(15.2062197472).margin(1e-3));
}

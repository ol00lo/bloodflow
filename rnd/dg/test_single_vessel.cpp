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
		return _points.back() / _points.size();
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
				ret.vals()[iaddr] -= v;
			}
		}

		return ret;
	}

	void save_vtk(const std::vector<double>& v, const std::string s) const {
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
	
		// Data
		fs << "POINT_DATA " << 2*n_elements() << std::endl;
		fs << "SCALARS data  double 1" << std::endl;
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
	double a_x0 = 0;
	double a_x1 = 0;
	double u_x0 = 0;
	double u_x1 = 0;
};

struct ProblemData{
	ProblemData(){
		beta = 4.0/3.0*sqrt(pi)*h*E/area0;
		amult = 4*sqrt(beta/2/rho);
		root4_a0 = sqrt(sqrt(area0));
	}

	static constexpr double pi = 3.1415926;

	// geometry parameters
	double L = 10;
	double area0 = pi*1e-4;

	// fluid parameters
	double rho = 1050;

	// vessel tissue parameters
	double h = 1.5*1e-3;
	double E = 400*1e3;

	// inflow conditions
	double q_inflow(double t) const{
		return 1e-6*exp(-1e4*(t-0.05)*(t-0.05));
	};

	// p(a)
	double pressure(double area) const{
		return 0 + beta*(std::sqrt(area) - std::sqrt(area0));
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
};

class IUpwindFluxCalculator{
public:
	virtual void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) = 0;
};

class InflowFluxCalculator: public IUpwindFluxCalculator{
public:
	InflowFluxCalculator(const FemGrid& grid, const ProblemData& data, std::function<double()> get_time, size_t cell):
		_data(data),
		_get_time(get_time),
		_cell_right(cell),
		_node_right(grid.tab_elem_nodes(cell)[0]),
		_sys(data.amult, data.root4_a0)
	{
		_a = data.area0;
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double t = _get_time();
		double q = _data.q_inflow(t);
		double w2_upw = _data.w2(area[_node_right], velocity[_node_right]);

		double area_upw = _a;
		double velo_upw = q / area_upw;
		_sys.set_qw(q, w2_upw);
		solve_nonlinear_system(_sys, area_upw, velo_upw);

		double fa = _data.flux_a(area_upw, velo_upw);
		double fu = _data.flux_u(area_upw, velo_upw);

		fluxes[_cell_right].a_x0 = fa;
		fluxes[_cell_right].u_x0 = fu;

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
				velo - _mult*(sqrt(sqrt(area)) - _root4_area0)
			};
		};
		std::array<double, 4> jac(double area, double velo) const override{
			return {
				velo, area,
				-0.25*_mult/area*area*area*sqrt(sqrt(area)), 1
			};
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

class MergingFluxCalculator: public IUpwindFluxCalculator{
public:
	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		_THROW_NOT_IMP_;
	}
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

TEST_CASE("Single vessel, constant properties", "[single-vessel-explicit]"){
	ProblemData data;
	double time = 0;

	FemGrid grid(data.L, 100, 1);
	double tau = grid.h()/2;
	std::vector<ElementBoundaryFluxes> fluxes(grid.n_elements());

	VtkUtils::TimeSeriesWriter writer("single-vessel");

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), 0.0);

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> upwind_flux_calculator(grid.n_points());
	upwind_flux_calculator[0].reset(new InflowFluxCalculator(grid, data, [&time](){ return time; }, 0));
	for (size_t i=1; i<grid.n_points()-1; ++i){
		upwind_flux_calculator[i].reset(new InternalFluxCalculator(grid, data, i-1, i));
	}
	upwind_flux_calculator.back().reset(new OutflowFluxCalculator(grid, data, grid.n_elements()-1));

	// prepare matrix solver
	CsrMatrix mass = grid.mass_matrix();
	CsrMatrix tran = grid.transport_matrix();
	AmgcMatrixSolver slv;
	slv.set_matrix(mass);
	
	while (time < 2.0){
		std::cout << time << std::endl;

		time += tau;
		
		// nodewise fluxes
		std::vector<double> flux_a(grid.n_nodes());
		std::vector<double> flux_u(grid.n_nodes());
		for (size_t i=0; i<grid.n_nodes(); ++i){
			flux_a[i] = data.flux_a(area[i], velocity[i]);
			flux_u[i] = data.flux_u(area[i], velocity[i]);
		}

		// upwind fluxes
		for (auto c: upwind_flux_calculator){
			c->compute(area, velocity, fluxes);
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
			rhs_a[node0] -= fluxes[ielem].a_x0;
			rhs_a[node1] += fluxes[ielem].a_x1;
			rhs_u[node0] -= fluxes[ielem].u_x0;
			rhs_u[node1] += fluxes[ielem].u_x1;
		}

		//solve
		slv.solve(rhs_a, area);
		slv.solve(rhs_u, velocity);

		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) grid.save_vtk(area, out_filename);
	}
}

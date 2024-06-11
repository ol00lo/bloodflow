#include "common_test.hpp"
#include "mat.hpp"
#include "vtk.hpp"
#include "nonlin_solver.hpp"

#define PI 3.1415926

namespace{

class FemGrid{
public:
	FemGrid(double length1, size_t n1,
	        double length2, size_t n2,
	        double length3, size_t n3,
	        size_t power)
	        : _power(power), _ne1(n1), _ne2(n2), _ne3(n3), _len1(length1), _len2(length2), _len3(length3){
	}

	double h(size_t ielem) const{
		if (ielem < _ne1){
			return _len1;
		} else if (ielem < _ne1 + _ne2){
			return _len2;
		} else {
			return _len3;
		}
	}

	size_t n_nodes() const {
		return n_elements()*(1 + _power);
	}

	size_t n_elements() const {
		return _ne1 + _ne2 + _ne3;
	}

	size_t n_points() const {
		return n_elements() + 1;
	}

	size_t n_local_bases() const {
		return _power + 1;
	}

	//double node(size_t i) const{
	//        return _nodes[i];
	//}

	double full_length() const{
		return _len1 + _len2 + _len3;
	}

	//size_t closest_node(double x) const {
	//        double t = 1e16;
	//        size_t ret = 0;
	//        for (size_t i=0; i<n_nodes(); ++i){
	//                double t1 = std::abs(_nodes[i] - x);
	//                if (t1 < t){
	//                        t = t1;
	//                        ret = i;
	//                }
	//        }
	//        return ret;
	//};

	CsrMatrix mass_matrix() const{
		CsrMatrix ret(stencil());
		std::vector<double> local = local_mass_matrix();

		for (size_t ielem=0; ielem<n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_nodes(ielem);

			for (size_t irow=0; irow<n_local_bases(); ++irow)
			for (size_t icol=0; icol<n_local_bases(); ++icol){
				double v = h(ielem)/2 * local[irow * n_local_bases() + icol];
				size_t iaddr = ret.get_address(lg[irow], lg[icol]);
				ret.vals()[iaddr] += v;
			}
		}

		return ret;
	}

	CsrMatrix block_transport_matrix() const{
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

	void save_vtk(const std::string s) const {
		std::ofstream fs(s);
		fs << "# vtk DataFile Version 3.0" << std::endl;
		fs << "DG" << std::endl;
		fs << "ASCII" << std::endl;
		fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		fs << "POINTS " << n_nodes() << " double" << std::endl;
		auto vtk_point_line = [&](double x0, double y0, double angle, double len, size_t n){
			double x1 = x0 + len * cos(angle);
			double y1 = y0 + len * sin(angle);
			double t;
			for (size_t i=0; i<n; ++i){
				t = (double)i/n;
				fs << (1-t)*x0 + t*x1 << " " << (1-t)*y0 + t*y1 << " 0" << std::endl;
				t = (double)(i+1)/n;
				fs << (1-t)*x0 + t*x1 << " " << (1-t)*y0 + t*y1 << " 0" << std::endl;
			}
		};
		vtk_point_line(0, 0, PI/2, _len1, _ne1);
		vtk_point_line(0, _len1, PI/3, _len2, _ne2);
		vtk_point_line(0, _len2, 2*PI/3, _len3, _ne3);

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
private:
	const size_t _power;
	const size_t _ne1, _ne2, _ne3;
	const double _len1, _len2, _len3;
	mutable CsrStencil _stencil;
	mutable std::vector<std::vector<size_t>> _tab_point_nodes;
	//std::vector<double> _points;
	//std::vector<double> _nodes;

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
		if (_tab_point_nodes.size() == 0){
			_tab_point_nodes.resize(n_points());
			for (size_t i=0; i<_ne1; ++i){
				size_t p1 = i;
				size_t p2 = i+1;
				_tab_point_nodes[p1].push_back(2*i);
				_tab_point_nodes[p2].push_back(2*i+1);
			}
			for (size_t i=_ne1; i<_ne2; ++i){
				size_t p1 = i;
				size_t p2 = i+1;
				_tab_point_nodes[p1].push_back(2*i);
				_tab_point_nodes[p2].push_back(2*i+1);
			}
			for (size_t i=_ne2; i<_ne3; ++i){
				size_t p1 = (i == _ne2) ? _ne1 : i;
				size_t p2 = i+1;
				_tab_point_nodes[p1].push_back(2*i);
				_tab_point_nodes[p2].push_back(2*i+1);
			}
		}
		return _tab_point_nodes[ipoint];
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


struct ProblemData{
	static constexpr double pi = PI;

	ProblemData(){
		recompute();
	}

	void recompute(){
		beta = 4.0/3.0*sqrt(pi)*h*E/area0;
		amult = 4*sqrt(beta/2/rho);
		root4_a0 = sqrt(sqrt(area0));
		visc_coef = -2*(profile_order+2)*mu*pi/rho;
	}

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

struct ElementBoundaryFluxes{
	double upwind_area_x0 = 0;
	double upwind_area_x1 = 0;
	double upwind_velo_x0 = 0;
	double upwind_velo_x1 = 0;
};

class IUpwindFluxCalculator{
public:
	virtual void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) = 0;
};

class InflowAUFluxCalculator: public IUpwindFluxCalculator{
public:
	InflowAUFluxCalculator(std::function<std::pair<double,double>()> getau, size_t cell):
		_cell_right(cell),
		_getau(getau)
	{
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
		double a, u;
		std::tie(a, u) = _getau();
		fluxes[_cell_right].upwind_area_x0 = a;
		fluxes[_cell_right].upwind_velo_x0 = u;
	}
private:
	const size_t _cell_right;
	std::function<std::pair<double,double>()> _getau;
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

class Bifurcation3FluxCalculator: public IUpwindFluxCalculator{
public:
	Bifurcation3FluxCalculator(
			const FemGrid& grid,
			const ProblemData& data_left,
			const ProblemData& data_right1,
			const ProblemData& data_right2,
			size_t elem_left,
			size_t elem_right1,
			size_t elem_right2,
			double eps=1e-12)
	: _elem1(elem_left), _elem2(elem_right1), _elem3(elem_right2),
	  _node1(grid.tab_elem_nodes(elem_left)[1]),
	  _node2(grid.tab_elem_nodes(elem_right1)[0]),
	  _node3(grid.tab_elem_nodes(elem_right2)[0]),
	  _data1(data_left), _data2(data_right1), _data3(data_right2),
	  _eps(eps), _sys(_data1, _data2, _data3)
	{
	}

	void compute(const std::vector<double>& area, const std::vector<double>& velocity, std::vector<ElementBoundaryFluxes>& fluxes) override{
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
	struct NonlinearSystem: public INonlinearSystem6{

		NonlinearSystem(const ProblemData& data1, const ProblemData& data2, const ProblemData& data3)
			:_data1(data1), _data2(data2), _data3(data3){}

		std::array<double, 6> f(double area1, double velo1, double area2, double velo2, double area3, double velo3) const override{
			return{
				_data1.flux_a(area1, velo1) - _data2.flux_a(area2, velo2) - _data3.flux_a(area3, velo3),
				_data1.flux_u(area1, velo1) - _data2.flux_u(area2, velo2),
				_data1.flux_u(area1, velo1) - _data2.flux_u(area3, velo3),
				_data1.w1(area1, velo1) - _w1,
				_data2.w2(area2, velo2) - _w2,
				_data3.w2(area3, velo3) - _w3
			};
		}
		std::array<double, 36> jac(double area1, double velo1, double area2, double velo2, double area3, double velo3) const override{
			return {
				velo1, area1, -velo2, -area2, -velo3, -area3,
				0.5*_data1.beta/_data1.rho/sqrt(area1), velo1, -0.5*_data2.beta/_data2.rho/sqrt(area2), -velo2, 0, 0,
				0.5*_data1.beta/_data1.rho/sqrt(area1), velo1, 0, 0, -0.5*_data3.beta/_data3.rho/sqrt(area3), -velo3,
				sqrt(_data1.beta/2/_data1.rho)*std::pow(area1, -0.75), 1, 0, 0, 0, 0,
				0, 0, -sqrt(_data2.beta/2/_data2.rho)*std::pow(area2, -0.75), 1, 0, 0,
				0, 0, 0, 0, -sqrt(_data3.beta/2/_data3.rho)*std::pow(area3, -0.75), 1
			};
		}

		void set_www(double w1, double w2, double w3){
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
	const ProblemData& _data1, _data2, _data3;
	const double _eps;
	NonlinearSystem _sys;

	double _a1, _a2, _a3;
	double _u1, _u2, _u3;
};


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
		_tran = _grid.block_transport_matrix();

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

}

TEST_CASE("Bifurcated vessel", "[bifurcated-vessel]"){
	ProblemData data1;
	data1.area0 = ProblemData::pi * 1e-4 / 4;
	data1.rho = 1000;
	data1.h = 1;
	data1.E = 3.0 * data1.area0 * 324.97 / 4.0 / std::sqrt(ProblemData::pi);
	data1.recompute();

	ProblemData data2;
	data2.area0 = ProblemData::pi * 1e-4 / 6.0 / 4.0;
	data2.rho = 1000;
	data2.h = 1;
	data2.E = 3.0 * data1.area0 * 796.02 / 4.0 / std::sqrt(ProblemData::pi);
	data2.recompute();


	double time = 0;
	double theta = 0.5;
	double L = 0.2;
	size_t k3 = 10;

	FemGrid grid(L, k3, L, k3, L, k3, 1);
	std::vector<int> cell_types(grid.n_elements(), 2);
	for (size_t i=0; i<k3; ++i){
		cell_types[i] = 1;
	}
	double tau = L/k3/50000;

	//size_t monitoring_node1 = grid.closest_node(0.25 * L);
	//size_t monitoring_node2 = grid.closest_node(0.5 * L);
	//size_t monitoring_node3 = grid.closest_node(0.75 * L);
	//std::vector<double> monitor1, monitor2, monitor3;

	std::vector<const ProblemData*> data(grid.n_elements());
	for (size_t icell=0; icell<grid.n_elements(); ++icell){
		data[icell] =  (cell_types[icell] == 1) ? &data1 : &data2;
	}

	std::vector<std::shared_ptr<IUpwindFluxCalculator>> flux_calculators(grid.n_points());
	// inflow
	auto getau = [&time]()->std::pair<double, double>{
		constexpr double time0 = 0.05;
		constexpr double C = 5000;
		constexpr double a = ProblemData::pi/4;

		double t = time - time0;
		double u = 0.01 * exp(-C*t*t);

		return {a, u};
	};
	flux_calculators[0].reset(new InflowAUFluxCalculator(getau, 0));
	// internal
	for (size_t i=1; i<k3; ++i) flux_calculators[i].reset(new InternalFluxCalculator(grid, data1, i-1, i));
	for (size_t i=k3+1; i<2*k3; ++i) flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i-1, i));
	for (size_t i=2*k3+1; i<3*k3; ++i) flux_calculators[i].reset(new InternalFluxCalculator(grid, data2, i-1, i));
	// outflow
	flux_calculators[2*k3].reset(new OutflowFluxCalculator(grid, data2, 2*k3-1));
	flux_calculators[3*k3].reset(new OutflowFluxCalculator(grid, data2, 3*k3-1));
	// bifurcation
	flux_calculators[k3].reset(new Bifurcation3FluxCalculator(grid, data1, data2, data2, k3-1, k3, 2*k3));

	// assembler
	Assembler assem(grid, data, flux_calculators);

	std::vector<double> velocity(grid.n_nodes(), 0.0);
	std::vector<double> area(grid.n_nodes(), data1.area0);

	AmgcMatrixSolver slv;

	VtkUtils::TimeSeriesWriter writer("bifurcated-vessel");
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
		std::cout << "TIME=" << time;
		std::cout << "  U=" << getau().second << std::endl;

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
		//monitor1.push_back(pressure[monitoring_node1]);
		//monitor2.push_back(pressure[monitoring_node2]);
		//monitor3.push_back(pressure[monitoring_node3]);
	}

	std::vector<double> pressure = assem.pressure(area);
	CHECK(pressure[20] == Approx(1521.4710994706).margin(1e-3));

	//time_value_vtk(tau, monitor1, "monitor1.vtk");
	//time_value_vtk(tau, monitor2, "monitor2.vtk");
	//time_value_vtk(tau, monitor3, "monitor3.vtk");
}

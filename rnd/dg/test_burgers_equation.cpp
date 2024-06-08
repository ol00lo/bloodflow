#include "common_test.hpp"
#include "mat.hpp"
#include <fstream>
#include "vtk.hpp"

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

	CsrMatrix mass_matrix() const{
		CsrMatrix ret(stencil());
		std::vector<double> local = local_mass_matrix();

		for (size_t ielem=0; ielem<n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_global_bases(ielem);

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
			std::vector<size_t> lg = tab_elem_global_bases(ielem);

			// block diagonal
			for (size_t irow=0; irow<n_local_bases(); ++irow)
			for (size_t icol=0; icol<n_local_bases(); ++icol){
				double v = local[irow * n_local_bases() + icol];
				size_t iaddr = ret.get_address(lg[irow], lg[icol]);
				ret.vals()[iaddr] -= v;
			}
			
			// upwind coupling
			if (ielem > 0){
				// left
				size_t iaddr = ret.get_address(lg[0], lg[0] - 1);
				ret.vals()[iaddr] -= 1;
			}
			{
				// right
				size_t iaddr = ret.get_address(lg[1], lg[1]);
				ret.vals()[iaddr] += 1;
			}
		}

		return ret;
	}

	CsrMatrix block_transport_matrix() const{
		CsrMatrix ret(stencil());
		std::vector<double> local = local_transport_matrix();

		for (size_t ielem=0; ielem<n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_global_bases(ielem);

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

	CsrMatrix coupled_transport_matrix() const{
		CsrMatrix ret(stencil());

		for (size_t ielem=0; ielem<n_elements(); ++ielem){
			std::vector<size_t> lg = tab_elem_global_bases(ielem);

			// upwind coupling
			if (ielem > 0){
				// left
				size_t iaddr = ret.get_address(lg[0], lg[0] - 1);
				ret.vals()[iaddr] -= 1;
			}
			{
				// right
				size_t iaddr = ret.get_address(lg[1], lg[1]);
				ret.vals()[iaddr] += 1;
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
			std::vector<size_t> lg = tab_elem_global_bases(ielem);
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
				std::vector<size_t> lg = tab_elem_global_bases(ielem);

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

	std::vector<size_t> tab_elem_global_bases(size_t ielem) const{
		std::vector<size_t> ret {2*ielem, 2*ielem+1};
		if (_power > 1){
			_THROW_NOT_IMP_;
		}
		return ret;
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

double exact(double x, double t=0){
	return x/(t+1);
}

}


TEST_CASE("Inviscid Burgers equation, explicit", "[Burgers-inviscid-explicit]"){
	FemGrid grid(1.0, 10, 1);
	double tau = grid.h()/2;

	CsrMatrix mass = grid.mass_matrix();
	CsrMatrix transport = grid.transport_matrix();
	CsrMatrix lhs = mass;
	// left boundary condition
	lhs.set_unit_row(0);

	AmgcMatrixSolver slv;
	slv.set_matrix(lhs);

	std::vector<double> u(grid.n_nodes());
	for (size_t i=0; i<grid.n_nodes(); ++i){
		u[i] = exact(grid.node(i));
	}
	VtkUtils::TimeSeriesWriter writer("burgers-transport");
	//writer.set_time_step(0.1);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) grid.save_vtk(u, out_filename);

	double time = 0;
	while (time + tau < 2.0 + 1e-6){
		time += tau;
		std::cout << time << std::endl;
		// assemble rhs
		std::vector<double> flux(grid.n_nodes());
		for (size_t i=0; i<grid.n_nodes(); ++i){
			flux[i] = u[i]*u[i]/2;
		}
		std::vector<double> fa = transport.mult_vec(flux);
		std::vector<double> rhs = mass.mult_vec(u);
		for (size_t i=0; i<rhs.size(); ++i){
			rhs[i] -= tau*fa[i];
		}

		// left boundary condition
		rhs[0] = 0.0;

		slv.solve(rhs, u);

		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) grid.save_vtk(u, out_filename);
	}
	CHECK(u.back() == Approx(0.3274573644).margin(1e-6));
}

TEST_CASE("Inviscid Burgers equation, implicit", "[Burgers-inviscid-implicit][amg]"){
	FemGrid grid(1.0, 10, 1);
	double tau = grid.h()/2;

	CsrMatrix mass = grid.mass_matrix();
	CsrMatrix transport = grid.transport_matrix();
	std::vector<double> load_vector = mass.mult_vec(std::vector<double>(grid.n_nodes(), 1.0));

	AmgcMatrixSolver slv(10'000, 1e-14);

	std::vector<double> u(grid.n_nodes());
	for (size_t i=0; i<grid.n_nodes(); ++i){
		u[i] = exact(grid.node(i));
	}
	VtkUtils::TimeSeriesWriter writer("burgers-transport");
	//writer.set_time_step(0.1);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) grid.save_vtk(u, out_filename);

	size_t maxit = 1000;
	double maxeps = 1e-14;
	double time = 0;
	while (time + tau < 2.0 + 1e-6){
		time += tau;
		std::cout << time << std::endl;

		// assemble rhs
		std::vector<double> rhs = mass.mult_vec(u);
		rhs[0] = 0.0;

		for (size_t it=0; it < maxit; ++it){
			// assemble lhs
			CsrMatrix lhs = mass;
			for (size_t i=0; i<lhs.n_rows(); ++i){
				for (size_t a=lhs.addr()[i]; a < lhs.addr()[i+1]; ++a){
					size_t col = lhs.cols()[a];
					lhs.vals()[a] += tau * u[col] / 2 * transport.vals()[a];
				}
			}
			// left boundary condition
			lhs.set_unit_row(0);

			// compute residual
			if (it > 0){
				std::vector<double> r = lhs.mult_vec(u);
				double err = 0;
				for (size_t i=0; i<r.size(); ++i){
					double diff = lhs.mult_vec(i, u) - rhs[i];
					err += diff*diff*load_vector[i];
				}
				err = std::sqrt(err);
				if (err < maxeps){
					std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				}
				if (it == maxit -1){
					std::cout << "Failed to converge with err=" << err << std::endl;
				}
			}

			slv.set_matrix(lhs);
			slv.solve(rhs, u);
		}

		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) grid.save_vtk(u, out_filename);
	}
	CHECK(u.back() == Approx(0.3396867099).margin(1e-6));
}

TEST_CASE("Inviscid Burgers equation, implicit, iterflux", "[Burgers-inviscid-implicit-iterflux][amg]"){
	FemGrid grid(1.0, 10, 1);
	double tau = grid.h()/5;

	CsrMatrix mass = grid.mass_matrix();
	CsrMatrix block_transport = grid.block_transport_matrix();
	CsrMatrix coupled_transport = grid.coupled_transport_matrix();
	std::vector<double> load_vector = mass.mult_vec(std::vector<double>(grid.n_nodes(), 1.0));

	AmgcMatrixSolver slv(10'000, 1e-14);

	std::vector<double> u(grid.n_nodes());
	for (size_t i=0; i<grid.n_nodes(); ++i){
		u[i] = exact(grid.node(i));
	}
	VtkUtils::TimeSeriesWriter writer("burgers-transport");
	//writer.set_time_step(0.1);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) grid.save_vtk(u, out_filename);

	size_t maxit = 10000;
	double maxeps = 1e-14;
	double time = 0;
	while (time + tau < 2.0 + 1e-6){
		time += tau;
		std::cout << time << std::endl;

		// assemble rhs
		std::vector<double> rhs1 = mass.mult_vec(u);

		for (size_t it=0; it < maxit; ++it){
			// assemble lhs
			CsrMatrix lhs = mass;
			for (size_t i=0; i<lhs.n_rows(); ++i){
				for (size_t a=lhs.addr()[i]; a < lhs.addr()[i+1]; ++a){
					size_t col = lhs.cols()[a];
					lhs.vals()[a] += tau * u[col] / 2 * block_transport.vals()[a];
				}
			}

			//assemble rhs
			std::vector<double> rhs(rhs1);
			std::vector<double> u2(grid.n_nodes());
			for (size_t i=0; i<u2.size(); ++i) u2[i] = u[i]*u[i];

			//// var1: goes to rhs
			//for (size_t i=0; i<grid.n_nodes(); ++i){
			//        rhs[i] -= tau * coupled_transport.mult_vec(i, u2) / 2.0;
			//}
	
			// var2: goes to the lhs diagonal
			for (size_t i=0; i<lhs.n_rows(); ++i){
				double denum = (u[i] == 0) ? 1e-16 : u[i];
				double v = tau * coupled_transport.mult_vec(i, u2) / 2.0 / denum;
				for (size_t a=lhs.addr()[i]; a < lhs.addr()[i+1]; ++a){
					size_t col = lhs.cols()[a];
					if (col == i){
						lhs.vals()[a] += v;
						break;
					}
				}
			}

			//// var3: goes to the both sides
			//for (size_t i=0; i<lhs.n_rows(); ++i){
			//        double v = coupled_transport.mult_vec(i, u2) / 2.0;
			//        size_t il, ir;
			//        if (i == 0 || i == lhs.n_rows() - 1){
			//                il = ir = i;
			//        } else if (i%2 == 0){
			//                il = i - 1;
			//                ir = i;
			//        } else {
			//                il = i;
			//                ir = i+1;
			//        }
			//        for (size_t a=lhs.addr()[i]; a < lhs.addr()[i+1]; ++a){
			//                size_t col = lhs.cols()[a];
			//                if (col == ir){
			//                        double d = (u[col] == 0) ? 1e-16 : u[col];
			//                        lhs.vals()[a] += tau * v / 2 / d;
			//                }
			//                if (col == il){
			//                        double d = (u[col] == 0) ? 1e-16 : u[col];
			//                        lhs.vals()[a] += tau * v / 2 / d;
			//                }
			//        }
			//}
	
			// left boundary condition
			lhs.set_unit_row(0);
			rhs[0] = 0;

			// compute residual
			if (it > 0){
				std::vector<double> r = lhs.mult_vec(u);
				double err = 0;
				for (size_t i=0; i<r.size(); ++i){
					double diff = lhs.mult_vec(i, u) - rhs[i];
					err += diff*diff*load_vector[i];
				}
				err = std::sqrt(err);
				if (err < maxeps){
					std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				}
				if (it == maxit -1){
					std::cout << "Failed to converge with err=" << err << std::endl;
				}
			}

			slv.set_matrix(lhs);
			slv.solve(rhs, u);
		}

		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) grid.save_vtk(u, out_filename);
	}
	CHECK(u.back() == Approx(0.3360909317).margin(1e-6));
}

TEST_CASE("Inviscid Burgers equation, cn", "[Burgers-inviscid-cn][amg]"){
	FemGrid grid(1.0, 10, 1);
	double tau = grid.h()/2;
	double theta = 0.5;

	CsrMatrix mass = grid.mass_matrix();
	CsrMatrix transport = grid.transport_matrix();
	std::vector<double> load_vector = mass.mult_vec(std::vector<double>(grid.n_nodes(), 1.0));

	AmgcMatrixSolver slv(10'000, 1e-14);

	std::vector<double> u(grid.n_nodes());
	for (size_t i=0; i<grid.n_nodes(); ++i){
		u[i] = exact(grid.node(i));
	}
	VtkUtils::TimeSeriesWriter writer("burgers-transport");
	//writer.set_time_step(0.1);
	std::string out_filename = writer.add(0);
	if (!out_filename.empty()) grid.save_vtk(u, out_filename);

	size_t maxit = 1000;
	double maxeps = 1e-14;
	double time = 0;
	while (time + tau < 2.0 + 1e-6){
		time += tau;
		std::cout << time << std::endl;

		// assemble rhs
		std::vector<double> rhs = mass.mult_vec(u);
		std::vector<double> flux(grid.n_nodes());
		for (size_t i=0; i<grid.n_nodes(); ++i){
			flux[i] = u[i]*u[i]/2;
		}
		flux = transport.mult_vec(flux);
		for (size_t i=0; i<grid.n_nodes(); ++i){
			rhs[i] -= (1-theta)*tau*flux[i];
		}
		rhs[0] = 0.0;

		for (size_t it=0; it < maxit; ++it){
			// assemble lhs
			CsrMatrix lhs = mass;
			for (size_t i=0; i<lhs.n_rows(); ++i){
				for (size_t a=lhs.addr()[i]; a < lhs.addr()[i+1]; ++a){
					size_t col = lhs.cols()[a];
					lhs.vals()[a] += theta * tau * u[col] / 2 * transport.vals()[a];
				}
			}
			// left boundary condition
			lhs.set_unit_row(0);

			// compute residual
			if (it > 0){
				std::vector<double> r = lhs.mult_vec(u);
				double err = 0;
				for (size_t i=0; i<r.size(); ++i){
					double diff = lhs.mult_vec(i, u) - rhs[i];
					err += diff*diff*load_vector[i];
				}
				err = std::sqrt(err);
				if (err < maxeps){
					std::cout << "converged in " << it << " iterations" << std::endl;
					break;
				}
				if (it == maxit -1){
					std::cout << "Failed to converge with err=" << err << std::endl;
				}
			}

			slv.set_matrix(lhs);
			slv.solve(rhs, u);
		}

		std::string out_filename = writer.add(time);
		if (!out_filename.empty()) grid.save_vtk(u, out_filename);
	}
	CHECK(u.back() == Approx(0.3335733261).margin(1e-6));
}

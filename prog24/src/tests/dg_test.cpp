#include "bflow/graph_grid.hpp"
#include "bflow/macros.hpp"
#include "bflow/matrix.hpp"
#include "bflow/matrix_solver.hpp"
#include "bflow/time_series_writer.hpp"
#include "bflow/vessel_graph.hpp"
#include "bflow/graph_grid.hpp"
#include "bflow/vtk.hpp"
#include "catch.hpp"
#include "bflow/graph_grid.hpp"
#include <fstream>
using namespace bflow;

class FemGrid
{
public:
    FemGrid(double length, size_t n, size_t power) : _power(power)
    {
        _points.push_back(0);
        _nodes.push_back(0);
        for (size_t i = 0; i < n; ++i)
        {
            _points.push_back(double(i + 1) / n * length);
            _nodes.push_back(_points.back());
            _nodes.push_back(_points.back());
        }
        _nodes.pop_back();
        if (power > 1)
        {
            _THROW_NOT_IMP_;
        }
    }
    FemGrid(const GraphGrid& grid) : _power(grid.n_midnodes)
    {
        if (_power != 1)
        {
            throw std::runtime_error("Only power=1 is allowed");
        }
        if (grid.n_cells() != 1)
        {
            throw std::runtime_error("Only single vessel grids are allowed");
        }
        _points.resize(grid.n_points(), 0);
        double x = 0;
        for (size_t i = 0; i < grid.n_elem(); ++i)
        {
            x += grid.find_cell_length(i);
            _points[i + 1] = x;
        }
        _nodes.resize(grid.n_nodes(), 0);
        for (size_t inode = 0; inode < _nodes.size(); ++inode)
        {
            size_t ipoint = (inode + 1) / 2;
            _nodes[inode] = _points[ipoint];
            _f_vec.push_back(grid.find_cell_length(0)/2);
        }
    }

    double h() const
    {
        return _points.back() / _points.size();
    }

    size_t n_nodes() const
    {
        return _nodes.size();
    }

    size_t n_elements() const
    {
        return _points.size() - 1;
    }

    size_t n_points() const
    {
        return _points.size();
    }

    size_t n_local_bases() const
    {
        return _power + 1;
    }

    double node(size_t i) const
    {
        return _nodes[i];
    }

    CsrMatrix mass_matrix() const
    {
        CsrMatrix ret(stencil());
        std::vector<double> local = local_mass_matrix();

        for (size_t ielem = 0; ielem < n_elements(); ++ielem)
        {
            std::vector<int> lg = tab_elem_global_bases(ielem);

            for (size_t irow = 0; irow < n_local_bases(); ++irow)
                for (size_t icol = 0; icol < n_local_bases(); ++icol)
                {
                    double v = h() / 2 * local[irow * n_local_bases() + icol];
                    size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                    ret.vals()[iaddr] += v;
                }
        }

        return ret;
    }

    CsrMatrix transport_matrix() const
    {
        CsrMatrix ret(stencil());
        std::vector<double> local = local_transport_matrix();

        for (size_t ielem = 0; ielem < n_elements(); ++ielem)
        {
            std::vector<int> lg = tab_elem_global_bases(ielem);

            // block diagonal
            for (size_t irow = 0; irow < n_local_bases(); ++irow)
                for (size_t icol = 0; icol < n_local_bases(); ++icol)
                {
                    double v = local[irow * n_local_bases() + icol];
                    size_t iaddr = ret.find_index(lg[irow], lg[icol]);
                    ret.vals()[iaddr] -= v;
                }

            // upwind coupling
            if (ielem > 0)
            {
                // left
                size_t iaddr = ret.find_index(lg[0], lg[0] - 1);
                ret.vals()[iaddr] -= 1;
            }
            {
                // right
                size_t iaddr = ret.find_index(lg[1], lg[1]);
                ret.vals()[iaddr] += 1;
            }
        }

        return ret;
    }

    std::vector<double> load_vector() const
    {
        return _f_vec;
    }

    void save_vtk(const std::vector<double>& v, const std::string s) const
    {
        std::ofstream fs(s);
        fs << "# vtk DataFile Version 3.0" << std::endl;
        fs << "DG" << std::endl;
        fs << "ASCII" << std::endl;
        fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
        fs << "POINTS " << n_nodes() << " double" << std::endl;
        for (const double& point : _nodes)
        {
            fs << point << " 0 0" << std::endl;
            fs << std::endl;
        }

        // Cells
        fs << "CELLS  " << n_elements() << "   " << 3 * n_elements() << std::endl;
        for (size_t ielem = 0; ielem < n_elements(); ++ielem)
        {
            std::vector<int> lg = tab_elem_global_bases(ielem);
            fs << 2 << " " << lg[0] << " " << lg[1] << std::endl;
        }
        fs << "CELL_TYPES  " << n_elements() << std::endl;
        for (size_t i = 0; i < n_elements(); ++i)
            fs << 3 << std::endl;

        // Data
        fs << "POINT_DATA " << 2 * n_elements() << std::endl;
        fs << "SCALARS data  double 1" << std::endl;
        fs << "LOOKUP_TABLE default" << std::endl;
        for (size_t i = 0; i < 2 * n_elements(); ++i)
        {
            fs << v[i] << std::endl;
        }
        fs.close();
    }

private:
    const int _power;
    mutable CsrMatrix _stencil;
    std::vector<double> _points;
    std::vector<double> _nodes;
    std::vector<double> _f_vec;

    CsrMatrix stencil() const
    {
        if (_stencil.n_rows() < 1)
        {
            std::vector<std::set<int>> lod(n_nodes());
            // block diagonal
            for (size_t ielem = 0; ielem < n_elements(); ++ielem)
            {
                std::vector<int> lg = tab_elem_global_bases(ielem);

                for (size_t irow = 0; irow < n_local_bases(); ++irow)
                    for (size_t icol = 0; icol < n_local_bases(); ++icol)
                    {
                        lod[lg[irow]].insert(lg[icol]);
                    }
            }
            // coupling
            for (int ipoint = 0; ipoint < n_points() - 1; ++ipoint)
            {
                std::vector<int> nodes = tab_point_nodes(ipoint);
                if (nodes.size() == 2)
                {
                    lod[nodes[0]].insert(nodes[1]);
                    lod[nodes[1]].insert(nodes[0]);
                }
            }
            _stencil.set_stencil(lod);
        }
        return _stencil;
    }

    std::vector<int> tab_point_nodes(int ipoint) const
    {
        if (ipoint == 0)
            return {0};
        else if (ipoint == n_points() - 1)
            return {2 * ipoint};
        else
            return {2 * ipoint - 1, 2 * ipoint};
    }

    std::vector<int> tab_elem_global_bases(int ielem) const
    {
        std::vector<int> ret{2 * ielem, 2 * ielem + 1};
        if (_power > 1)
        {
            _THROW_NOT_IMP_;
        }
        return ret;
    }

    std::vector<double> local_mass_matrix() const
    {
        if (_power == 1)
        {
            return {2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0};
        }
        else
        {
            _THROW_NOT_IMP_;
        }
    }

    std::vector<double> local_transport_matrix() const
    {
        if (_power == 1)
        {
            return {-0.5, -0.5, 0.5, 0.5};
        }
        else
        {
            _THROW_NOT_IMP_;
        }
    }
};

namespace
{
double exact(double x, double t = 0)
{
    constexpr double eps = 1e-6;
    if (t > 0)
        return exact(x - t);
    else
        return (x > eps && x < 0.8 - eps) ? 1.0 : 0.0;
}

double norm2(FemGrid grid, std::vector<double> u, double time)
{   
    auto lv = grid.load_vector();
    double I = 0;
    double gamma = 0;
    for (size_t i = 0; i < grid.n_nodes() ; ++i)
    {
        double r = u[i] - exact(grid.node(i), time);
        I +=  lv[i] * (r * r);
        gamma += grid.h();
    }
    return std::sqrt(I / gamma);
}
}

TEST_CASE("Transport equation, upwind", "[upwind-transport]")
{
    // Legacy test. Can be removed
    FemGrid grid(3.0, 30, 1);  
    double tau = grid.h() / 2;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.transport_matrix();
    CsrMatrix lhs = mass;
    CsrMatrix rhs_mat = mass;
    for (size_t i = 0; i < mass.n_nonzeros(); ++i)
    {
        lhs.vals()[i] += tau / 2.0 * transport.vals()[i];
        rhs_mat.vals()[i] -= tau / 2 * transport.vals()[i];
    }
    // left boundary condition
    lhs.set_unit_row(0);

    // matrix solver
    AmgcMatrixSolver slv;
    slv.set_matrix(lhs);

    // initial conditions
    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact(grid.node(i));
    }
    TimeSeriesWriter writer("upwind-transport");
    writer.set_time_step(0.1);
    std::string out_filename = writer.add(0);
    if (!out_filename.empty())
        grid.save_vtk(u, out_filename);


    std::vector<double> L(grid.n_points(),0.5);

    double time = 0;
    while (time < 2.0)
    {
        std::cout << time << std::endl;

        // assemble rhs
        std::vector<double> rhs = rhs_mat.mult_vec(u);
        // left boundary condition
        rhs[0] = 0.0;

        slv.solve(rhs, u);
        time += tau;

        std::string out_filename = writer.add(time);
        if (!out_filename.empty())
            grid.save_vtk(u, out_filename);
    }
    CHECK(u[50] == Approx(1.071884924).margin(1e-6));
}

TEST_CASE("Transport equation, upwind2", "[upwind-transport2]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {3.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1);
    FemGrid grid(grid1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    double tau = grid.h() / 2;

    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.transport_matrix();
    CsrMatrix lhs = mass;
    CsrMatrix rhs_mat = mass;
    for (size_t i = 0; i < mass.n_nonzeros(); ++i)
    {
        lhs.vals()[i] += tau / 2.0 * transport.vals()[i];
        rhs_mat.vals()[i] -= tau / 2 * transport.vals()[i];
    }
    // left boundary condition
    lhs.set_unit_row(0);

    // matrix solver
    AmgcMatrixSolver slv;
    slv.set_matrix(lhs);

    // initial conditions
    std::vector<double> u(grid.n_nodes());
    for (size_t i = 0; i < grid.n_nodes(); ++i)
    {
        u[i] = exact(grid.node(i));
    }
    TimeSeriesWriter writer("upwind-transport1");
    //writer.set_time_step(0.1);
    std::string out_filename = writer.add(0);
    if (!out_filename.empty())
        grid.save_vtk(u, out_filename);


    double time = 0;
    while (time < 2.0)
    {
        // assemble rhs
        std::vector<double> rhs = rhs_mat.mult_vec(u);
        // left boundary condition
        rhs[0] = 0.0;

        slv.solve(rhs, u);
        time += tau;
        
        std::cout << time << "  ";
        std::cout << norm2(grid, u, time) << std::endl;

        std::string out_filename = writer.add(time);
        if (!out_filename.empty())
            grid.save_vtk(u, out_filename);
    }
    CHECK(u[50] == Approx(1.071884924).margin(1e-6));
    CHECK(norm2(grid, u, time) == Approx(0.102352).margin(1e-4));
}

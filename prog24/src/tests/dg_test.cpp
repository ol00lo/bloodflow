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
#include <sstream>
#include <iomanip>
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
        _h = length / n;
        if (power > 1)
        {
            _THROW_NOT_IMP_;
        }
    }
    FemGrid(const GraphGrid& grid, std::vector<Point2> nodes_coo) : _power(grid._power)
    {
        if (_power > 4)
        {
            throw std::runtime_error("Only power=1,2,3,4 is allowed");
        }
        if (grid.n_edges() != 1)
        {
            throw std::runtime_error("Only single vessel grids are allowed");
        }
        _points.resize(grid.n_points(), 0.0);
        double x = 0;
        for (size_t i = 0; i < grid.n_elem(); ++i)
        {
            x += grid.find_cell_length(i);
            _points[i + 1] = x;
        }
        _nodes.resize(grid.n_nodes()-(_power-1)*(_points.size()-1), 0);
        _h = grid.find_cell_length(0);

        for (size_t inode = 0; inode < _nodes.size(); ++inode)
        {
            size_t ipoint = (inode + 1) / 2;
            _nodes[inode] = _points[ipoint];
            _f_vec.push_back(grid.find_cell_length(0) / 2);
        }
        double l = _h / _power;
        for (size_t i = 0; i < _points.size()-1; i++)
        {
            double p = _points[i] + l;
            for (size_t j = 1; j < _power; j++)
            {
                _nodes.push_back(p);
                p += l;
            }
        }
        
    }

    double h() const
    {
        return _h;
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

private:
    const int _power;
    mutable CsrMatrix _stencil;
    std::vector<double> _points;
    std::vector<double> _nodes;
    std::vector<double> _f_vec;
    double _h;

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
        std::vector<int> ret;
        if (_power == 1)
        {
            ret = {2 * ielem, 2 * ielem + 1};
        }
        else if (_power == 2)
        {
            ret = {2 * ielem, 2 * ielem + 1, 2 * int(n_elements()) + ielem};
        }
        else if (_power == 3)
        {
            ret = {
                2 * ielem, 
                2 * ielem + 1,
                2 * int(n_elements()) + 2 * ielem,
                2 * int(n_elements()) + 2 * ielem + 1 };
        }
        else if (_power == 4)
        {
            ret = {
                2 * ielem,
                2 * ielem + 1,
                2 * int(n_elements()) + 3 * ielem,
                2 * int(n_elements()) + 3 * ielem + 1,
                2 * int(n_elements()) + 3 * ielem + 2 };
        }
        else
        {
            _THROW_NOT_IMP_;
        }
        return ret;
    }

    std::vector<double> local_mass_matrix() const
    {
        if (_power == 1)
        {
            return {2.0 / 3.0, 1.0 / 3.0,
                    1.0 / 3.0, 2.0 / 3.0};
        }
        else if (_power == 2)
        {
            return {4.0 / 15.0, -1.0 / 15.0,  2.0 / 15.0,
                   -1.0 / 15.0,  4.0 / 15.0,  2.0 / 15.0,
                    2.0 / 15.0,  2.0 / 15.0, 16.0 / 15.0};
        }
        else if (_power == 3)
        {
            return {16.0 / 105.0,  19.0 / 840.0,  33.0 / 280.0, -3.0 / 70,
                    19.0 / 840.0,  16.0 / 105.0,  -3.0 / 70.0,  33.0 / 280.0,
                    33.0 / 280.0 , -3.0 / 70.0,   27.0 / 35.0, -27.0 / 280.0,
                    -3.0 / 70.0,   33.0 / 280.0, -27.0 / 280.0, 27.0 / 35.0};
        }
        else if (_power == 4)
        {
            return {292.0/2835, -29.0/2835, 296.0/2835,  -58.0/945,    8.0/405,
                    -29.0/2835, 292.0/2835,   8.0/405,   -58.0/945,  296.0/2835,
                    296.0/2835,   8.0/405,  256.0/405,  -128.0/945,  256.0/2835,
                    -58.0/945,  -58.0/945, -128.0/945,   208.0/315, -128.0/945,
                      8.0/405,  296.0/2835, 256.0/2835, -128.0/945,  256.0/405};
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
        else if (_power == 2)
        {
            return {-1.0 / 2.0, 1.0 / 6.0, -2.0 / 3.0,
                    -1.0 / 6.0, 1.0 / 2.0, 2.0 / 3.0, 
                    2.0 / 3.0, -2.0 / 3.0, 0.0};
        }
        else if (_power == 3)
        {
            return {-1.0 / 2.0,   -7.0 / 80.0, -57.0 / 80.0, 3.0 / 10.0,
                     7.0 / 80.0,   1.0 / 2.0,   -3.0 / 10.0, 57.0 / 80.0,
                    57.0 / 80.0,   3.0 / 10.0,     0.0,     -81.0 / 80.0,
                    -3.0 / 10.0, -57.0 / 80.0,  81.0 / 80.0,    0.0};
        }
        else if (_power == 4)
        {
            return { -1.0/2,    107.0/1890, -736.0/945, 134.0/315, -64.0/315,
                   -107.0/1890,   1.0/2,      64.0/315,-134.0/315, 736.0/945,
                   736.0/945,  -64.0/315,      0.0,   -352.0/315, 512.0/945,
                   -134.0/315,  134.0/315,   352.0/315,    0.0,   -352.0/315,
                     64.0/315, -736.0/945,  -512.0/945, 352.0/315,    0.0 };        
        }
        else
        {
            _THROW_NOT_IMP_;
        }
    }
};

namespace{
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
    for (size_t i = 0; i < grid.n_points() ; ++i)
    {
        double r = u[i] - exact(grid.node(i), time);
        I +=  lv[i] * (r * r);
        gamma += grid.h();
    }
    return std::sqrt(I / gamma);
}


void print_matrix_full(const CsrMatrix& mat, std::ostream& s = std::cout)
{
    for (size_t i = 0; i < mat.n_rows(); i++)
    {
        for (size_t j = 0; j < mat.n_rows(); j++)
        {
            if (mat.is_in_stencil(i, j) == false)
            {
                s << std::setw(6) << "*";
            }
            else
            {
                s << std::setw(6) << round(mat.value(i, j) * 1000) / 1000;
            }
        }
        s << std::endl;
    }
}
}
TEST_CASE("Transport equation, upwind", "[upwind-transport]")
{
    // Legacy test. Can be removed
    FemGrid grid(3.0, 30, 1);  
    double tau = grid.h() / 3;

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
    TimeSeriesWriter writer("upwind-transport-1");
    writer.set_time_step(0.1);
    std::string out_filename = writer.add(0);
    if (!out_filename.empty())
        //grid.save_vtk(u, out_filename);


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
        //if (!out_filename.empty())
            //grid.save_vtk(u, out_filename);
    }
    CHECK(u[50] == Approx(1.071884924).margin(1e-6));
}


TEST_CASE("Transport equation, upwind2", "[upwind-transport2]")
{
    std::vector<std::vector<int>> node = {{0}, {0}};
    std::vector<double> ed = {4.0};
    VesselGraph gr1(node, ed);
    GraphGrid grid1(gr1, 0.1, 1);
    std::vector<Point2> nodes_coo = generate_nodes_coo(gr1);
    FemGrid grid(grid1, nodes_coo);
    double tau = grid.h() / 2;
    CsrMatrix mass = grid.mass_matrix();
    CsrMatrix transport = grid.transport_matrix();
    CsrMatrix lhs = mass;
    CsrMatrix rhs_mat = mass;
    for (size_t i = 0; i < mass.n_nonzeros(); ++i)
    {
        //ab
        //lhs.vals()[i] += 3.0*tau / 2.0 * transport.vals()[i];
        //rhs_mat.vals()[i] += tau / 2 * transport.vals()[i];
        //kn
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

    NonstatGridSaver saver(grid1, nodes_coo, "dg1kn");
    saver.new_time_step(0);
    saver.save_vtk_point_data(u, "data");
    for (double t = tau; t <= 2.0; t += tau)
    {
        std::vector<double> rhs = rhs_mat.mult_vec(u);

        // left boundary condition
        rhs[0] =0.0;

        slv.solve(rhs, u);

        std::cout << t << "  ";
        std::cout << norm2(grid, u, t) << std::endl;

        saver.new_time_step(t);
        saver.save_vtk_point_data(u, "data");
    }
    //CHECK(u[50] == Approx(1.071884924).margin(1e-6));
    CHECK(norm2(grid, u, 2.0) == Approx(0.102352).margin(1e-4));
}

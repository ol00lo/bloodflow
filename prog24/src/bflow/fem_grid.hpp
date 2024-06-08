#ifndef FEM_GRID_HPP
#define FEM_GRID_HPP
#include "bflow/graph_grid.hpp"
namespace bflow
{
class FemGrid
{
public:
    FemGrid(const GraphGrid& grid, std::vector<Point2> nodes_coo);
    double h() const;
    size_t n_nodes() const;
    size_t n_elements() const;
    size_t n_points() const;
    size_t n_local_bases() const;
    double node(size_t i) const;
    std::vector<double> load_vector() const;
    std::vector<int> tab_elem_nodes(int ielem) const;
    double full_length() const;
    CsrMatrix mass_matrix() const;
    CsrMatrix transport_matrix() const;
    CsrMatrix block_transport_matrix() const;
    CsrMatrix coupled_transport_matrix() const;

private:
    const int _power;
    mutable CsrMatrix _stencil;
    std::vector<double> _points;
    std::vector<double> _nodes;
    std::vector<double> _f_vec;
    double _h;
    double _len = 0;
    CsrMatrix stencil() const;
    std::vector<double> local_mass_matrix() const;
    std::vector<double> local_transport_matrix() const;
    std::vector<int> tab_point_nodes(int ipoint) const;
    void fill_f_vec();
};
} // namespace bflow
#endif
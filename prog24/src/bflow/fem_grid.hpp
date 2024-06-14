#ifndef FEM_GRID_HPP
#define FEM_GRID_HPP
#include "bflow/graph_grid.hpp"
namespace bflow
{
class FemGrid
{
public:
    FemGrid(const GraphGrid& grid);
    double h(int ielem) const;
    size_t n_nodes() const;
    size_t n_elements() const;
    size_t n_points() const;
    size_t n_local_bases() const;
    double node(size_t i) const;
    std::vector<double> load_vector() const;
    std::vector<int> tab_elem_nodes(int ielem) const;
    double full_length() const;
    // integral[ phi_j phi_i dx]
    CsrMatrix mass_matrix() const;
    // integral[ (phi_j) (d phi_i / dx) dx]
    CsrMatrix block_transport_matrix() const;
    // integral[ u (phi_j) (d phi_i / dx) dx]
    CsrMatrix block_u_transport_matrix(const std::vector<double>& u) const;
    // phi_j phi_i (x1) - phi_j phi_i (x0)
    CsrMatrix coupling_transport_matrix() const;
    size_t closest_node(double x) const;
    size_t closest_node(size_t section, double p) const;
    size_t tab_node_elem(size_t inode) const;
private:
    const int _power;
    mutable CsrMatrix _stencil;
    std::vector<double> _points;
    std::vector<double> _nodes;
    std::vector<double> _f_vec;
    std::vector<double> _h;
    double _len = 0;
    CsrMatrix stencil() const;
    std::vector<double> local_mass_matrix() const;
    std::vector<double> local_u_transport_matrix() const;
    std::vector<double> local_transport_matrix() const;
    std::vector<int> tab_point_nodes(int ipoint) const;
    void fill_f_vec();
    std::vector<std::array<int, 2>> _nodes_by_edge;
};
} // namespace bflow
#endif

#ifndef VESSEL_GRAPH_HPP
#define VESSEL_GRAPH_HPP

#include <array>
#include <iostream>
#include <map>
#include <vector>

namespace bflow
{
class VesselGraph
{
public:
    VesselGraph(const std::vector<std::vector<int>>& tab_node_edges, const std::vector<double>& edge_lens);
    std::array<int, 2> tab_edge_node(int iedge) const;
    std::vector<int> tab_node_edge(int inode) const;
    int find_edge_by_nodes(int inode, int jnode) const;
    double find_length(int iedge) const;
    int n_nodes() const;
    int n_edges() const;

private:
    std::vector<std::vector<int>> _nodes_edges;
    std::vector<double> _edge_lengths;
    int _n_edges;
};
} // namespace bflow

#endif
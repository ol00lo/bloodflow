#include <iostream>
#include <vector>
#include <map>
#include <array>

namespace bflow {
	class VesselGraph {
	public:
		VesselGraph(const std::vector<std::vector<int>>& node, const std::vector<double>& len);
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

	class GraphGrid {
	public:
		GraphGrid(const VesselGraph& _graph, double h);
		int n_points() const;
		int n_cells() const;
		std::vector<int> tab_point_cells(int point) const;
		std::array<int, 2> tab_cell_points(int cell) const;
		int find_edge_by_cell(int cell) const;
		std::vector<int> points_by_edge(int edge) const;
	private:
		std::vector<std::vector<int>> _points;
		std::vector<std::vector<double>> _cells;
	};
}

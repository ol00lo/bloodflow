#include "vessel_graph.hpp"
#include <algorithm>


bflow::VesselGraph::VesselGraph(const std::vector<std::vector<int>>& edges, const std::vector<double>& len) {
	if (edges.size() < 2 || len.size() < 1 || edges.size()-len.size() != 1) throw std::runtime_error("Incorrect input"); 
	
	_nodes_edges = edges;
	_edge_lengths = len;

	int i = 0;
	std::vector<int> all_edges;
	for (size_t node = 0; node < _nodes_edges.size(); ++node) {
		for (int edge = 0; edge < _nodes_edges[node].size(); ++edge) {
			all_edges.push_back(_nodes_edges[node][edge]);
			i++;
		}
	}

	_n_edges = i / 2;

	if (len.size() != _n_edges) throw std::runtime_error("Incorrect input");

	for (size_t edge = 0; edge < all_edges.size(); ++edge) {
		int cnt = std::count(cbegin(all_edges), cend(all_edges), all_edges[edge]);
		if (cnt != 2) throw std::runtime_error("Incorrect input");
	}
}

std::array<int, 2> bflow::VesselGraph::tab_edge_node(int iedge) const {
	std::array<int, 2> nodes = { -1, -1 };

	int i = -1;
	for (size_t node = 0; node < _nodes_edges.size(); ++node) 
	for (size_t edge = 0; edge < _nodes_edges[node].size(); ++edge) {
			if (_nodes_edges[node][edge] == iedge) {
				i++;
				if (i > 1) throw std::runtime_error("Incorrect input");
				nodes[i] = node;
			}	
	}
	if (nodes[0] == -1 || nodes[1] == -1) {
		throw std::runtime_error("Incorrect input");
	}
	return nodes;
}

std::vector<int> bflow::VesselGraph::tab_node_edge(int inode) const {
	if (inode >= _nodes_edges.size()) {
		throw std::runtime_error("Out of graph");
	}
	return _nodes_edges[inode];
}

int bflow::VesselGraph::find_edge_by_nodes(int inode, int jnode) const {
	if (inode >= _nodes_edges.size()|| jnode >= _nodes_edges.size()) {
		throw std::runtime_error("Out of graph");
	}
	std::vector<int> iedges=_nodes_edges[inode];
	std::vector<int> jedges=_nodes_edges[jnode];
	for(size_t i = 0; i < iedges.size(); ++i)
	for (size_t j = 0; j < jedges.size(); ++i) {
			if (jedges[j] == iedges[i]) return iedges[i];
	}
	return -1;
}

double bflow::VesselGraph::find_length(int iedge) const {
	if (iedge >= _edge_lengths.size()) {
		throw std::runtime_error("Out of graph");
	}
	return _edge_lengths[iedge];
}

int bflow::VesselGraph::n_nodes() const {
	return _nodes_edges.size();
}

int bflow::VesselGraph::n_edges() const {
	return _n_edges;
}
#include "vessel_graph.hpp"
#include <algorithm>
#include <set>

using namespace bflow;

// VESSEL GRAPH

VesselGraph::VesselGraph(const std::vector<std::vector<int>>& edges, const std::vector<double>& len) {
	if (edges.size() < 2 || len.size() < 1 || edges.size() - len.size() != 1) throw std::runtime_error("Incorrect input");

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

std::array<int, 2> VesselGraph::tab_edge_node(int iedge) const {
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

std::vector<int> VesselGraph::tab_node_edge(int inode) const {
	if (inode >= _nodes_edges.size()) {
		throw std::runtime_error("Out of graph");
	}
	return _nodes_edges[inode];
}

int VesselGraph::find_edge_by_nodes(int inode, int jnode) const {
	if (inode >= _nodes_edges.size() || jnode >= _nodes_edges.size()) {
		throw std::runtime_error("Out of graph");
	}
	std::vector<int> iedges = _nodes_edges[inode];
	std::vector<int> jedges = _nodes_edges[jnode];
	for (size_t i = 0; i < iedges.size(); ++i)
		for (size_t j = 0; j < jedges.size(); ++i) {
			if (jedges[j] == iedges[i]) return iedges[i];
		}
	return -1;
}

double VesselGraph::find_length(int iedge) const {
	if (iedge >= _edge_lengths.size()) {
		throw std::runtime_error("Out of graph");
	}
	return _edge_lengths[iedge];
}

int VesselGraph::n_nodes() const {
	return _nodes_edges.size();
}

int VesselGraph::n_edges() const {
	return _n_edges;
}


// GRAPH GRID

GraphGrid::GraphGrid(const VesselGraph& graph, double h) {
	_points.resize(graph.n_edges());
	_cells.resize(graph.n_edges());
	int n_nodes = graph.n_nodes();

	for (int edge = 0; edge < graph.n_edges(); ++edge) {

		std::array<int, 2> boundary = graph.tab_edge_node(edge);
		double len = graph.find_length(edge);
		int n_cells = round(len / h);
		double h_cell = 1.0 * len / n_cells;

		_points[edge].push_back(boundary[0]);
		for (size_t cells = 1; cells < n_cells; ++cells) {
			_points[edge].push_back(n_nodes);
			_cells[edge].push_back(h_cell * cells);
			n_nodes++;
		}
		_points[edge].push_back(boundary[1]);
		_cells[edge].push_back(len);
	}
}

int GraphGrid::n_points() const {
	std::set<int> n_points;
	for (size_t i = 0; i < _points.size(); ++i)
		for (size_t j = 0; j < _points[i].size(); ++j) {
			n_points.insert(_points[i][j]);
		}
	return n_points.size();
}

int GraphGrid::n_cells() const {
	int n_cells = 0;
	for (size_t i = 0; i < _cells.size(); ++i)
		for (size_t j = 0; j < _cells[i].size(); ++j) {
			n_cells++;
		}
	return n_cells;
}

std::array<int, 2> GraphGrid::tab_cell_points(int cell) const {
	std::array<int, 2> points = { -1,-1 };
	int i_cell = 0;

	for (size_t i = 0; i < _points.size(); ++i)
		for (size_t j = 1; j < _points[i].size(); ++j) {
			if (i_cell == cell) {
				points[0] = _points[i][j - 1];
				points[1] = _points[i][j];
			}
			i_cell++;
		}
	return points;
}

std::vector<int> GraphGrid::tab_point_cells(int point) const {
	std::vector<int> cells;
	int i_cell = 0;
	for (size_t i = 0; i < _points.size(); ++i) {
		for (size_t j = 0; j < _points[i].size(); ++j) {
			if (_points[i][j] == point) {
				if (j == _points[i].size() - 1) {
					cells.push_back(i_cell - 1);
				}
				if (j == 0) {
					cells.push_back(i_cell);
				}
				if (j != 0 && j != _points[i].size() - 1) {
					cells.push_back(i_cell - 1);
					cells.push_back(i_cell);
				}
			}
			i_cell++;
		}
		i_cell--;
	}
	return cells;
}

int GraphGrid::find_edge_by_cell(int cell) const {
	int i_cell = 0;
	for (size_t i = 0; i < _points.size(); ++i) {
		for (size_t j = 0; j < _points[i].size(); ++j) {
			if (i_cell == cell) {
				return i;
			}
			i_cell++;
		}
		i_cell--;
	}
	return -1;
}

std::vector<int> GraphGrid::points_by_edge(int edge) const {
	return _points[edge];
}
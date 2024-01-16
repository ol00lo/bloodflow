#include "catch.hpp"
#include <vector>
#include <map>

class Graph {
public:
	Graph(std::vector<std::vector<int>> node, std::vector<double> len);
	//void add_node(std::vector<int> node, std::vector<double> len);
	std::map<int, std::vector<int>> edge_node();
	int find_edge(int a, int b);
	double find_length(int a, int b);
	double find_lenght(int a);
	int n_nodes();
	int n_edges();
private:
	std::map<int, std::vector<int>> _nodes;
	std::map<int, double> _edges;
};

Graph::Graph(std::vector<std::vector<int>> node, std::vector<double> len) {
	for (size_t nod = 0; nod < node.size(); ++nod) {
		_nodes.insert( {nod, node[nod]} );
	}
	for (size_t le = 0; le < len.size(); ++le) {
		_edges.insert({ le, len[le] });
	}
}

std::map<int, std::vector<int>> Graph::edge_node() {
	std::map<int, std::vector<int>> e_n;
	int i = 0;
	for (auto node = _nodes.begin(); node != _nodes.end(); ++node) {
		for (int j = 0; j < node->second.size(); ++j) {
			std::vector<int> nodes;
			if (node->second[j] > node->first) {
				nodes.push_back(node->first);
				nodes.push_back(node->second[j]);
				e_n.insert({ i, nodes });
				i++;
			}
		}
		
	}
	return e_n;
}

int Graph::find_edge(int a, int b) {
	std::vector<int> nodes1 = { a, b };
	std::vector<int> nodes2 = { b, a };

	std::map<int, std::vector<int>> e_n = edge_node();
	for (auto ed = e_n.begin(); ed != e_n.end(); ++ed) {
		if (ed->second == nodes1 || ed->second == nodes2) {
			return ed->first;
		}
	}
	return -1;
}

double Graph::find_lenght(int a) {
	auto c = _edges.find(a);
	return c->second;
}

double Graph::find_length(int a, int b) {
	int c = find_edge(a, b);
	return find_lenght(c);
}

int Graph::n_nodes() {
	return _nodes.size();
}

int Graph::n_edges() {
	int i = 0;
	for (auto it = _nodes.begin(); it != _nodes.end(); ++it) {
		for (int j = 0; j < it->second.size(); ++j) {
			i++;
		}
	}
	return i / 2;
}

TEST_CASE("basic graph functionality", "[graph]") {
	std::vector<std::vector<int>> node = { {1,2,3}, {0, 4, 5}, {0, 6, 7}, {0}, {1}, {1, 8, 9}, {2, 10, 11, 12}, {2}, {5}, {5}, {6}, {6}, {6} };
	std::vector<double> ed = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
	Graph gr1(node, ed);
	CHECK(gr1.n_edges() == 12);
	CHECK(gr1.edge_node().size() == gr1.n_edges());
	CHECK(gr1.find_edge(7,2) == 6);
	CHECK(gr1.find_lenght(7) == 8.0);
	CHECK(gr1.find_length(7, 2) == 7.0);
}
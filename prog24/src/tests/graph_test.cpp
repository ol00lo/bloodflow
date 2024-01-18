#include "catch.hpp"
#include <vector>
#include <map>
#include "vessel_graph.hpp"


TEST_CASE("basic graph functionality", "[graph]") {
	std::vector<std::vector<int>> node = { {0,1,2}, {0}, {1, 3, 4, 5}, {2, 6, 7}, {3}, {4}, {5}, {6, 8, 9, 10}, {7}, {8}, {9}, {10} };
	std::vector<double> ed = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 };
	VesselGraph gr1(node, ed);
	CHECK(gr1.n_edges() == ed.size());
	CHECK(gr1.find_length(7) == 8.0);
	CHECK(gr1.tab_edge_node(4)[0] == 2);
	CHECK(gr1.tab_edge_node(4)[1] == 5);
}

TEST_CASE("improper graph definition", "[improper_graph]") {
	CHECK_THROWS(VesselGraph({ {0,1,2} }, { 1.0, 2.0 }));
	CHECK_THROWS(VesselGraph({ {0,1,2}, {0}, {1}, {2} }, { 1.0, 2.0 }));
	CHECK_THROWS(VesselGraph({ {0,1,2}, {0}, {1}, {1} }, { 1.0, 2.0, 3.0 }));
	CHECK_THROWS(VesselGraph({ {0,1,2}, {0}, {1} }, { 1.0, 2.0 }));
	CHECK_THROWS(VesselGraph({ {0,1,2,2}, {0}, {1} }, { 1.0, 2.0 }));
	CHECK_THROWS(VesselGraph({ {0,1,2,2}, {0}, {1} }, { 1.0, 2.0, 3.0 }));
}
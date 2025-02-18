
#include <iostream>

#include "branch_n_bound_seq.hpp"
#include "branching_strategy.hpp"

#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "dsatur_color.hpp"

#include "csr_graph.hpp"
#include "dimacs.hpp"
#include "test_common.hpp"

int main() {
	Dimacs dimacs;
	std::string file_name = "10_vertices_graph.col";

	if (!dimacs.load(file_name.c_str())) {
		std::cout << dimacs.getError() << std::endl;
		return 1;
	}

	CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);

	RandomBranchingStrategy branching_strategy(graph->GetNumVertices());
	FastCliqueStrategy clique_strategy;
	DSaturColorStrategy color_strategy;

	BranchNBoundSeq solver(branching_strategy, clique_strategy,
			       color_strategy, "log.txt");

	int chromatic_number = solver.Solve(*graph, 10000, 100000);

	std::cout << "Chromatic number: " << chromatic_number << std::endl;

	return 0;
}


#include <iostream>

#include "branch_n_bound_seq.hpp"
#include "branching_strategy.hpp"

#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "dsatur_color.hpp"
#include "advanced_color.hpp"

#include "csr_graph.hpp"
#include "dimacs.hpp"
#include "test_common.hpp"

int main(int argc, char* argv[]) {
	Dimacs dimacs;

	int cutoff_time 	  = 60;
	int cutoff_iterations = 100000;

	std::string file_name = "queen10_10.col";
	if ( argc >= 2 ) {
		file_name = argv[1];

		if ( argc >= 3 ) {
			cutoff_time = std::stoi(argv[2]);

			if ( argc >= 4 ) {
				cutoff_iterations = std::stoi(argv[3]);
			}
		} 
	}

	if (!dimacs.load(file_name.c_str())) {
		std::cout << dimacs.getError() << std::endl;
		return 1;
	}

	CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);

	NeighboursBranchingStrategy branching_strategy;

	FastCliqueStrategy clique_strategy;

	DSaturColorStrategy base_color_strategy;
	GreedySwapRecolorStrategy recolor_strategy;
	ColorNRecolorStrategy first_color_strategy(base_color_strategy, recolor_strategy);
	InactiveColorStrategy second_color_strategy;
	InterleavedColorStrategy color_strategy(first_color_strategy, second_color_strategy, 
										 	1, 10);

	BranchNBoundSeq solver(branching_strategy, clique_strategy,
			       color_strategy, "log.txt");

	int chromatic_number = solver.Solve(*graph, cutoff_time, cutoff_iterations);

	std::cout << "Chromatic number: " << chromatic_number << std::endl;

	return 0;
}

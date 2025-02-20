
#include <iostream>

#include "branch_n_bound_par.hpp"
#include "branching_strategy.hpp"

#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "dsatur_color.hpp"

#include "csr_graph.hpp"
#include "dimacs.hpp"
#include "test_common.hpp"

#include <mpi.h>

//TODO: case p=1 cause problems
int main(int argc, char** argv) {
	Dimacs dimacs;
	std::string file_name = "10_vertices_graph.col";

	if (!dimacs.load(file_name.c_str())) {
		std::cout << dimacs.getError() << std::endl;
		return 1;
	}

	CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);


	NeighboursBranchingStrategy branching_strategy;
	FastCliqueStrategy clique_strategy;
	GreedyColorStrategy color_strategy;

	BranchNBoundPar solver(branching_strategy, clique_strategy,
			       color_strategy, "log_master.txt", "log_branches.txt");

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	if (provided < MPI_THREAD_MULTIPLE) {
		std::cerr << "MPI does not support full multithreading!" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	int chromatic_number = solver.Solve(*graph, 5, 100000);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if (my_rank == 0)
		std::cout << "Chromatic number: " << chromatic_number << std::endl;

	MPI_Finalize();
	return 0;
}

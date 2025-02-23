
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

int main(int argc, char** argv) {
	Dimacs dimacs;
	std::string file_name = "10_vertices_graph.col";

	if (!dimacs.load(file_name.c_str())) {
		std::cout << dimacs.getError() << std::endl;
		return 1;
	}
	std::cout << "Succesfully read Graph." << std::endl;

	CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);

	NeighboursBranchingStrategy branching_strategy;
	FastCliqueStrategy clique_strategy;
	GreedyColorStrategy color_strategy;

	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	if (provided < MPI_THREAD_MULTIPLE) {
		std::cerr << "MPI does not support full multithreading!" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	double optimum_time;
	int chromatic_number = solver.Solve(*graph, optimum_time, 10);

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	BranchNBoundPar solver(branching_strategy, clique_strategy, color_strategy, "log_" + std::to_string(my_rank) + ".txt");
	
	if (my_rank == 0) {
		if (optimum_time == -1)
			std::cout << "Rank 0: Finalizing with timeout" << std::endl;
		else 
			std::cout << "Rank 0: Finalizing with time " << optimum_time << std::endl;
		std::cout << "Chromatic number: " << chromatic_number << std::endl;
	}
	MPI_Finalize();
	return 0;
}

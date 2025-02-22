
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
	CSRGraph* graph;
	NeighboursBranchingStrategy branching_strategy;
	FastCliqueStrategy clique_strategy;
	DSaturColorStrategy color_strategy;

	// Initialize MPI with multithreading enabled
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided < MPI_THREAD_MULTIPLE) {
		std::cerr << "MPI does not support full multithreading!" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	double optimum_time;
	// Test the solver multiple times on the same graph (since its not deterministic)
	int N_trials = 20;
	for (int i=0; i<N_trials; i++) {
		// Root process reads the graph
		if (my_rank == 0) {
			
			if (!dimacs.load(file_name.c_str())) {
				std::cout << dimacs.getError() << std::endl;
				return 1;
			}
			
			std::cout << "Succesfully read Graph." << std::endl;

			graph = CSRGraph::LoadFromDimacs(file_name);

		}
	
		BranchNBoundPar solver(branching_strategy, clique_strategy,
			color_strategy, "log_par");
		if (my_rank == 0) 
			std::cout << "Starting trial " << i << "..." << std::endl;
		int chromatic_number = solver.Solve(*graph, optimum_time, 15);
		if (my_rank == 0) {
			if (optimum_time == -1)
				std::cout << "Rank 0: Finalizing with timeout" << std::endl;
			else 
				std::cout << "Rank 0: Finalizing with time " << optimum_time << std::endl;
			std::cout << "Chromatic number: " << chromatic_number << std::endl;
		}
	}
	MPI_Finalize();
	return 0;
}

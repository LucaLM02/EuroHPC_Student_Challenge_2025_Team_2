
#include <iostream>
#include <fstream>

#include "balanced_branch_n_bound_par.hpp"
#include "branching_strategy.hpp"

#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "recolor.hpp"
#include "advanced_color.hpp"
#include "dsatur_color.hpp"

#include "csr_graph.hpp"
#include "dimacs.hpp"
#include "test_common.hpp"

#include <mpi.h>

#include <cstdlib> // For std::stoi

int main(int argc, char** argv) {
    // Default values
    int timeout = 60;
    int N_trials = 1;
    std::string folder_name = "graphs_instances/";
    std::string file_name;
    int expected_chi;

    // Check for required arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <file_name> <expected_chi> [timeout] [N_trials]\n";
        return 1;
    }

    file_name = folder_name + argv[1]; // Prefix folder to filename
    expected_chi = std::stoi(argv[2]);

    // Parse optional timeout argument
    if (argc > 3) {
        try {
            timeout = std::stoi(argv[3]); // Convert timeout to integer
            if (timeout <= 0) {
                std::cerr << "Error: Timeout must be a positive integer.\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid timeout argument.\n";
            return 1;
        }
    }

    // Parse optional N_trials argument
    if (argc > 4) {
        try {
            N_trials = std::stoi(argv[4]); // Convert N_trials to integer
            if (N_trials <= 0) {
                std::cerr << "Error: N_trials must be a positive integer.\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid N_trials argument.\n";
            return 1;
        }
    }


    Dimacs dimacs;
    CSRGraph* graph;
    NeighboursBranchingStrategy branching_strategy;
    FastCliqueStrategy clique_strategy;
    // GreedyColorStrategy greedy_color_strategy;
    // DSaturColorStrategy base_color_strategy;
    // GreedySwapRecolorStrategy recolor_strategy;
    // ColorNRecolorStrategy advanced_color_strategy(base_color_strategy, recolor_strategy);

    // InterleavedColorStrategy color_strategy(greedy_color_strategy, advanced_color_strategy, 5, 2);
    GreedyColorStrategy color_strategy;

	// Initialize MPI with multithreading enabled
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided < MPI_THREAD_MULTIPLE) {
		std::cerr << "MPI does not support full multithreading!" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	int my_rank, n_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

    if ( my_rank == 0 ) {
        std::cout << "Reading file: " << file_name << "\n";
        std::cout << "Using timeout: " << timeout << " seconds\n";
        std::cout << "Number of trials: " << N_trials << "\n";
    }

	double optimum_time;
	// Test the solver multiple times on the same graph (since its not deterministic)
	for (int i=0; i<N_trials; i++) {
		// Root process reads the graph
        if (!dimacs.load(file_name.c_str())) {
            std::cout << dimacs.getError() << std::endl;
            return 1;
        }

        graph = CSRGraph::LoadFromDimacs(file_name);

		if (my_rank == 0) {
			
            std::cout << "Succesfully read Graph." << std::endl;
		}
	
        BalancedBranchNBoundPar solver(branching_strategy, clique_strategy, color_strategy, "log" + std::to_string(my_rank) + ".txt");

		if (my_rank == 0) {
		std::cout << "Starting trial " << i << "..." << std::endl;
		}
		// Start the timer.
		auto start_time = MPI_Wtime();
		int chromatic_number = solver.Solve(*graph, optimum_time, timeout, expected_chi);
		auto end_time = MPI_Wtime();
		auto time = end_time - start_time;
		if (my_rank == 0) {
			std::cout << "Trial " << i << " took " << time << " seconds." << std::endl;
			if (optimum_time == -1)
				std::cout << "It was a timeout." << std::endl;
			else 
				std::cout << "Finalized prematurely with self-measured " << optimum_time << " seconds. " << std::endl;
			std::cout << "Chromatic number: " << chromatic_number << std::endl;

            {
                std::ofstream out("chromatic_number.txt");
                out << chromatic_number << std::endl;
                out.close();
            }

            std::ofstream out("requested_output.txt");
            out << "problem_instance_file_name "    << file_name << std::endl;
            out << "cmd line "                      << std::endl;
            out << "solver version "                << std::endl;
            out << "number_of_vertices "            << graph->GetNumVertices() << std::endl;
            out << "number_of_edges: "              << graph->GetNumEdges() << std::endl;
            out << "time_limit_sec "                << timeout << std::endl;
            out << "number_of_worker_processes "    << n_proc << std::endl;
            out << "number_of_cores_per_worker "    << 4 << std::endl;
            if ( optimum_time == -1 ) {
                out << "wall_time_sec "             << "> 10000" << std::endl;
                out << "is_within_time_limit "      << false << std::endl;
            } else {
                out << "wall_time_sec "             << optimum_time << std::endl;
                out << "is_within_time_limit "      << true << std::endl;
            }

            std::vector<unsigned short> colors = graph->GetFullColoring();
            unsigned max_color = 0;
            for ( int color : colors ) {
                if ( color > max_color ) {
                    max_color = color;
                }
            }

            TestFunctions::CheckColoring(*graph);
            std::cout << graph->HasEdge(6, 136) << std::endl;
            std::cout << graph->GetColor(6) << std::endl;
            std::cout << graph->GetColor(136) << std::endl;

            out << "number_of_colors "              << max_color << std::endl;
            std::vector<int> vertices = graph->GetVertices();
            for ( int vertex : vertices ) {
                out << vertex << " " << colors[vertex] << std::endl;
            }

		}
	}

	MPI_Finalize();
	return 0;
}

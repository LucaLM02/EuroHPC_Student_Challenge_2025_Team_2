#include <iostream>
#include <fstream>

#include "branch_n_bound_par.hpp"
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
    std::string output_file_name;
    unsigned short expected_chi;

    // Check for required arguments
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <file_name> <expected_chi> <timeout> <output file> <color_strategy_flag> <b&b_algorithm_flag>\n";
        return 1;
    }

    file_name = folder_name + argv[1]; // Prefix folder to filename
    expected_chi = static_cast<unsigned short>(std::stoi(argv[2]));

    // Parse optional timeout argument
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

    output_file_name = argv[4];

    int strategy_flag = std::stoi(argv[5]);
    ColorStrategy* color_strategy;
    GreedyColorStrategy greedy_color_strategy;
    DSaturColorStrategy base_color_strategy;
    GreedySwapRecolorStrategy recolor_strategy;
    ColorNRecolorStrategy advanced_color_strategy(base_color_strategy, recolor_strategy);
    InterleavedColorStrategy heavy_color_strategy(greedy_color_strategy, advanced_color_strategy, 5, 2);
    if ( strategy_flag == 0) {
        // light color strategy
        color_strategy = &greedy_color_strategy;
    } else if ( strategy_flag == 1 ) {
        // heavy color strategy
        color_strategy = &heavy_color_strategy;
    }

    Dimacs dimacs;
    CSRGraph* graph;
    NeighboursBranchingStrategy branching_strategy;
    FastCliqueStrategy clique_strategy;



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
    if (!dimacs.load(file_name.c_str())) {
        std::cout << dimacs.getError() << std::endl;
        return 1;
    }

    graph = CSRGraph::LoadFromDimacs(file_name);

    if (my_rank == 0) {
        
        std::cout << "Succesfully read Graph." << std::endl;
    }

    int b_n_b_flag = std::stoi(argv[6]);
    int chromatic_number;
    double start_time, end_time;
    if ( b_n_b_flag == 0) {
        BalancedBranchNBoundPar solver(branching_strategy, clique_strategy, *color_strategy, "log" + std::to_string(my_rank) + ".txt");
        start_time = MPI_Wtime();
        chromatic_number = solver.Solve(*graph, optimum_time, timeout, 10 ,expected_chi);
        end_time = MPI_Wtime();
    } else if ( b_n_b_flag == 1 ) {
        BranchNBoundPar solver(branching_strategy, clique_strategy, *color_strategy, "log" + std::to_string(my_rank) + ".txt");
        start_time = MPI_Wtime();
        chromatic_number = solver.Solve(*graph, optimum_time, timeout, 10 ,expected_chi);
        end_time = MPI_Wtime();
    }


    // Start the timer.
    auto time = end_time - start_time;
    if (my_rank == 0) {
        if (optimum_time == -1)
            std::cout << "It was a timeout." << std::endl;
        else 
            std::cout << "Finalized prematurely with self-measured " << optimum_time << " seconds. " << std::endl;
        std::cout << "Chromatic number: " << chromatic_number << std::endl;

        std::ofstream out(output_file_name);
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

        out << "number_of_colors "              << max_color << std::endl;
        std::vector<int> vertices = graph->GetVertices();
        for ( int vertex : vertices ) {
            out << vertex << " " << colors[vertex] << std::endl;
        }

    }

	MPI_Finalize();
	return 0;
}

#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cstdlib> // For std::stoi
#include <unordered_map>
#include <filesystem>

#include "branch_n_bound_par.hpp"
#include "branching_strategy.hpp"
#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "dsatur_color.hpp"
#include "csr_graph.hpp"
#include "dimacs.hpp"

int main(int argc, char** argv) {
    // Default values
    int timeout = 60;
    int sol_gather_period = 10;
    std::string file_name;

    // Check for required arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file_name> [timeout]\n";
        return 1;
    }

    file_name = argv[1]; // Get filename from arguments

    // Parse optional arguments
    if (argc > 2) {
        try {
            timeout = std::stoi(argv[2]); // Convert timeout to integer
            if (timeout <= 0) {
                std::cerr << "Error: Timeout must be a positive integer.\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid timeout argument.\n";
            return 1;
        }
        try {
            sol_gather_period = std::stoi(argv[3]); // Convert timeout to integer
            if (sol_gather_period <= 0) {
                std::cerr << "Error: Solution gathering period must be a positive integer.\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid sol_gather_period argument.\n";
            return 1;
        }
    }

    std::cout << "Reading file: " << file_name << "\n";
    std::cout << "Using timeout: " << timeout << " seconds\n";
    std::cout << "Using sol_gather_period: " << sol_gather_period << " seconds\n";

    // Load expected results from text file
    std::ifstream txt_file("expected_chi.txt");
    if (!txt_file.is_open()) {
        std::cerr << "Error: Could not open expected results text file.\n";
        return 1;
    }

    std::unordered_map<std::string, int> expected_results;
    std::string key;
    int value;
    while (txt_file >> key >> value) {
        expected_results[key] = value;
    }
    txt_file.close();

    // Check if the expected result for the given file is available
    std::string file_key = std::filesystem::path(file_name).filename().string();
    if (expected_results.find(file_key) == expected_results.end()) {
        std::cerr << "Error: No expected result found for the given file.\n";
        return 1;
    }

    int expected_chromatic_number = expected_results[file_key];

    Dimacs dimacs;
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

    std::string full_file_name = "graphs_instances/" + file_name; 
    // All processes read the graph, since they all start with it.
    if (!dimacs.load(full_file_name.c_str())) {
        std::cout << dimacs.getError() << std::endl;
        return 1;
    }

    std::cout << "Rank " << my_rank << ": Successfully read Graph " << file_name << std::endl;

    graph = CSRGraph::LoadFromDimacs(full_file_name);

    BranchNBoundPar solver(branching_strategy, clique_strategy, color_strategy, "logs/log_" + std::to_string(my_rank) + ".txt");

    // Start the timer.
    auto start_time = MPI_Wtime();
    double optimum_time;
    int chromatic_number = solver.Solve(*graph, optimum_time, timeout-0.05, sol_gather_period);
    auto end_time = MPI_Wtime();
    auto time = end_time - start_time;

    if (my_rank == 0) {
        std::cout << "Execution took " << time << " seconds." << std::endl;
        if (optimum_time == -1)
            std::cout << "It was a timeout." << std::endl;
        else
            std::cout << "Solve() finished prematurely measuring " << optimum_time << " seconds. " << std::endl;
        

        // Compare with expected chromatic number
        if (chromatic_number != expected_chromatic_number) {
            std::cout << "Failed: expected " << expected_chromatic_number << " but got " << chromatic_number << std::endl;
        } else {
            std::cout << "Chromatic number: " << chromatic_number << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}

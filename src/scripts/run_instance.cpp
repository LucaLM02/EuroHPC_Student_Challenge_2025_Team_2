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

/**
 * @brief Main function to run the graph coloring solver using the branch and bound method.
 *
 * This function initializes MPI, reads the graph instance from a file, and runs the BranchNBoundPar solver.
 * It handles command-line arguments for the input file, timeout, and solution gathering period.
 * The function compares the computed chromatic number with the expected result and prints the outcome and
 * computation time.
 */
int main(int argc, char** argv) {
    // Default values
    int timeout = 60;
    int sol_gather_period = 10;
    int balanced = 1;
    std::string file_name;

    // Check for required arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file_name> [--timeout=<timeout>] [--sol_gather_period=<period>] [--balanced=<0|1>]\n";
        return 1;
    }

    file_name = argv[1]; // Get filename from arguments

    // Parse optional arguments
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        std::istringstream iss(arg);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            try {
                if (key == "--timeout") {
                    timeout = std::stoi(value);
                    if (timeout <= 0) {
                        std::cerr << "Error: Timeout must be a positive integer.\n";
                        return 1;
                    }
                } else if (key == "--sol_gather_period") {
                    sol_gather_period = std::stoi(value);
                    if (sol_gather_period <= 0) {
                        std::cerr << "Error: Solution gathering period must be a positive integer.\n";
                        return 1;
                    }
                } else if (key == "--balanced") {
                    balanced = std::stoi(value);
                } else {
                    std::cerr << "Error: Unknown argument " << arg << "\n";
                    return 1;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error: Invalid value for argument " << key << ".\n";
                return 1;
            }
        } else {
            std::cerr << "Error: Invalid argument format " << arg << ".\n";
            return 1;
        }
    }

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

    // Output arguments
    if (my_rank == 0) {
        std::cout << "Reading file: " << file_name << "\n";
        std::cout << "Using timeout: " << timeout << " seconds\n";
        std::cout << "Using sol_gather_period: " << sol_gather_period << " seconds\n";
        std::cout << "Using balanced approach: " << balanced << "\n";
    }

    // Read the Graph
    std::string full_file_name = "graphs_instances/" + file_name; 
    // All processes read the graph, since they all start with it.
    if (!dimacs.load(full_file_name.c_str())) {
        std::cout << dimacs.getError() << std::endl;
        return 1;
    }
    graph = CSRGraph::LoadFromDimacs(full_file_name);
    std::cout << "Rank " << my_rank << ": Successfully read Graph " << file_name << std::endl;

    BranchNBoundPar solver(branching_strategy, clique_strategy, color_strategy, "logs/log_" + std::to_string(my_rank) + ".txt");
    BalancedBranchNBoundPar balanced_solver(branching_strategy, clique_strategy, color_strategy, "logs/log_" + std::to_string(my_rank) + ".txt");


    // Start the timer.
    auto start_time = MPI_Wtime();
    // Run.
    double optimum_time;    
    int chromatic_number;
    if (balanced) {
        chromatic_number = balanced_solver.Solve(*graph, optimum_time, timeout-0.05, sol_gather_period,  expected_chromatic_number);
    } else {
        chromatic_number = solver.Solve(*graph, optimum_time, timeout-0.05, sol_gather_period, expected_chromatic_number);
    }
    // Stop the timer.
    auto end_time = MPI_Wtime();
    auto time = end_time - start_time;

    // Output results
    if (my_rank == 0) {
        std::cout << "Execution took " << time << " seconds." << std::endl;
        if (optimum_time == -1)
            std::cout << "It was a timeout." << std::endl;
        else
            std::cout << "Solve() finished prematurely measuring " << optimum_time << " seconds. " << std::endl;
        
        // Compare with expected chromatic number
        if (chromatic_number != expected_chromatic_number) 
            std::cout << "Failed: expected " << expected_chromatic_number << " but got " << chromatic_number << std::endl;
        else 
            std::cout << "Suceeded: Chromatic number: " << chromatic_number << std::endl;
        
    }

    MPI_Finalize();
    return 0;
}

#ifndef BRANCH_N_BOUND_SEQ_HPP
#define BRANCH_N_BOUND_SEQ_HPP

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>

#include "bound_strategies.hpp"
#include "branching_strategy.hpp"
#include "graph.hpp"

/**
 * @brief Sequential branch-and-bound solver.
 *
 * This class implements a sequential branch-and-bound algorithm to solve for
 * the chromatic number of a graph.
 */
class BranchNBoundSeq {
       private:
	BranchingStrategy& _branching_strat;
	CliqueStrategy& _clique_strat;
	ColorStrategy& _color_strat;
	std::ofstream _log_file;

	/**
	 * @brief Checks if the solver has exceeded the timeout.
	 *
	 * @param start_time The start time of the solver.
	 * @param timeout_seconds The timeout duration in seconds.
	 * @return bool True if the timeout has been reached, false otherwise.
	 */
	static bool CheckTimeout(
	    const std::chrono::steady_clock::time_point& start_time,
	    int timeout_seconds);

	/**
	 * @brief Constructs a BranchNBoundSeq solver.
	 *
	 * @param timeout_seconds The maximum time (in seconds) the solver can
	 * run before it times out.
	 * @param iteration_threshold The maximum number of iterations without
	 * improvement before stopping.
	 */
       public:
	BranchNBoundSeq(BranchingStrategy& branching_strat,
			CliqueStrategy& clique_strat,
			ColorStrategy& color_strat,
			const std::string& log_file_path)
	    : _branching_strat(branching_strat),
	      _clique_strat(clique_strat),
	      _color_strat(color_strat) {
		_log_file.open(log_file_path);
		if (!_log_file.is_open()) {
			throw std::runtime_error("Failed to open log file: " +
						 log_file_path);
		}
	}

	/**
	 * @brief Solves a graph using the sequential branch-and-bound
	 * algorithm.
	 *
	 * This method computes and returns the chromatic number for the input
	 * graph.
	 *
	 * @param g The graph to solve.
	 */
	int Solve(Graph& g, int timeout_seconds = 60,
		  int iteration_threshold = 1000);

	// Destructor to close the log file
	~BranchNBoundSeq() {
		if (_log_file.is_open()) {
			_log_file << "Closing log file.\n";
			_log_file.close();
		}
	}
};

#endif	// BRANCH_N_BOUND_SEQ_HPP

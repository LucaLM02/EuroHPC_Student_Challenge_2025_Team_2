#ifndef BRANCH_N_BOUND_PAR_HPP
#define BRANCH_N_BOUND_PAR_HPP

#include <mpi.h>

#include <chrono>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <queue>
#include <utility>
#include <omp.h>
#include <mutex>
#include <atomic>

#include "common.hpp"
#include "graph.hpp"
#include "branching_strategy.hpp"
#include "clique_strategy.hpp"
#include "color.hpp"

class BranchNBoundPar {
       private:
	BranchingStrategy& _branching_strat;
	CliqueStrategy& _clique_strat;
	ColorStrategy& _color_strat;
	std::ofstream _log_file;

	/**
	 * @brief Logs a message to the log file.
	 *
	 * @param message The message to log.
	 */
	 void Log(const std::string& message, int depth, bool is_branching);

	/**
	 * @brief Logs a message to the log file.
	 *
	 * @param message The message to log.
	 */
	 void Log_par(const std::string& message, int depth, bool is_branching);

	/**
	 * @brief Creates a task for the OpenMP parallel region to execute.
	 *
	 */
	 void create_task(std::atomic<int>& active_tasks, Graph* current_G, int u, int v,
			 CliqueStrategy& _clique_strat, ColorStrategy& _color_strat,
			 std::vector<Branch>& new_branches, int task_type,
			 std::atomic<unsigned short>const &best_ub, int const &depth);

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

       public:
	/**
	 * @brief Constructs a BranchNBoundPar solver.
	 *
	 * @param timeout_seconds The maximum time (in seconds) the solver can
	 * run before it times out.
	 * @param iteration_threshold The maximum number of iterations without
	 * improvement before stopping.
	 */
	BranchNBoundPar(BranchingStrategy& branching_strat,
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

	// TODO: Either input vector of Graphs to solve or one at a time popped
	// from a Graph queue?
	int Solve(Graph& g, int timeout_seconds = 60,
		  int iteration_threshold = 1000);
};

#endif	// BRANCH_N_BOUND_PAR_HPP

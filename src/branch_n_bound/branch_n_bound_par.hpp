#ifndef BRANCH_N_BOUND_PAR_HPP
#define BRANCH_N_BOUND_PAR_HPP

#include <mpi.h>
#include <omp.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <utility>
#include <vector>
#include <thread>

#include "branching_strategy.hpp"
#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "common.hpp"
#include "graph.hpp"

class BranchNBoundPar {
       private:
	BranchingStrategy& _branching_strat;
	CliqueStrategy& _clique_strat;
	ColorStrategy& _color_strat;
	std::ofstream _log_file;

	void create_task(
		std::unique_ptr<Graph> current_G, int u, int v,
		CliqueStrategy& _clique_strat, ColorStrategy& _color_strat,
		std::atomic<unsigned short> & best_ub, int const& depth, int my_rank);

	/**
	 * @brief Logs a message to the log file in a openmp parallel section.
	 *
	 * @param message The message to log.
	 */
	void Log_par(const std::string& message, int depth, bool is_branching, int rank, int thread_id);

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

	int Solve(Graph& g, double &optimum_time, int timeout_seconds = 60);
};

#endif	// BRANCH_N_BOUND_PAR_HPP

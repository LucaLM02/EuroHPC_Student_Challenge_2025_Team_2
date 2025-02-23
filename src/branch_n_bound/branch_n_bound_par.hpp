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

using BranchQueue = std::priority_queue<Branch, std::vector<Branch>>;

class BranchNBoundPar {
	private:
		BranchingStrategy& _branching_strat;
		CliqueStrategy& _clique_strat;
		ColorStrategy& _color_strat;
		std::ofstream _log_file;
	
		// void create_task(
		// 	std::unique_ptr<Graph> current_G, int u, int v,
		// 	CliqueStrategy& _clique_strat, ColorStrategy& _color_strat,
		// 	std::atomic<unsigned short> & best_ub, int const& depth, int my_rank);
	
		/**
		 * @brief Logs a message to the log file in a openmp parallel section.
		 *
		 * @param message The message to log.
		 */
		void Log_par(const std::string& message, int depth);
	
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
		 * @brief Listens for termination signals (solution found or timeout) and terminates the execution of the process accordingly.
		 *
		 * @param my_rank The rank of the current process.
		 * @param p The total number of processes in the MPI communicator.
		 * @param global_start_time The global start time, used to check for timeouts.
		 * @param timeout_seconds The timeout duration (in seconds) after which the timeout signal is sent.
		 * @param optimum_time The time at which the optimum solution was found.
		 */
		void thread_0_terminator(int my_rank, int p, int global_start_time, int timeout_seconds, double &optimum_time);
	
		/**
		 * @brief Updates (gathers) best_ub from time to time.
		 *
		 * @param p The total number of processes in the MPI communicator.
		 * @param best_ub The best upper bound found so far.
		 * @param my_rank The rank of the current process.
		 */
		void thread_1_solution_gatherer(int p, std::atomic<unsigned short> &best_ub);
	
		/**
		 * @brief Employer thread employs workers by answering their work requests.
		 *
		 * @param queue_mutex Mutex to protect concurrent access to the work queue.
		 * @param queue The local work queue containing branches to be processed.
		 */
		void thread_2_employer(std::mutex& queue_mutex, BranchQueue& queue);
	
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
						const std::string& log_file_path_prefix)
			: _branching_strat(branching_strat),
			  _clique_strat(clique_strat),
			  _color_strat(color_strat) {
				int my_rank;
				MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
				_log_file.open(log_file_path_prefix + std::to_string(my_rank) + ".log");
				if (!_log_file.is_open()) {
					throw std::runtime_error("Failed to open log file: " +
					log_file_path_prefix + std::to_string(my_rank) + ".log");
			}
		}
	
		int Solve(Graph& g, double &optimum_time, int timeout_seconds = 60);
	};

#endif	// BRANCH_N_BOUND_PAR_HPP

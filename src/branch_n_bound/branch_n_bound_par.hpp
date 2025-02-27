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

		std::atomic<unsigned short> _best_ub = USHRT_MAX;
		std::mutex _best_branch_mutex;
		Branch _current_best;

		void ColorInitialGraph(Graph& initial_graph, const Branch& optimal_branch);

		/**
		 * @brief updates the _current_best branch to store the best graph with the best coloring
		 * 
		 * @param depth depth at which the graph was colored
		 * @param lb lower bound of the best colored graph
		 * @param ub upper bound of the best colored graph
		 * @param graph the best colored graph
		 */
		void UpdateCurrentBest(int depth, int lb, unsigned short ub, GraphPtr graph);

		/**
		 * @brief Logs a message to the log file in a MPI environment and OpenMP parallel section.
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
		void thread_0_terminator(int my_rank, int p, int global_start_time, int timeout_seconds, 
								 double &optimum_time, Graph& graph_to_color);
	
		/**
         * @brief Updates (gathers) best_ub from time to time.
         *
         * @param p The total number of processes in the MPI communicator.
         * @param best_ub The best upper bound found so far.
         * @param my_rank The rank of the current process.
         * @param sol_gather_period The period (in seconds) at which the best upper bound is gathered.
         */
		void thread_1_solution_gatherer(int p, std::atomic<unsigned short> &best_ub, int sol_gather_period);
	
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
         * @param branching_strat The branching strategy to use.
         * @param clique_strat The clique strategy to use.
         * @param color_strat The color strategy to use.
         * @param log_file_path The path to the log file.
         */
		 BranchNBoundPar(BranchingStrategy& branching_strat,
			CliqueStrategy& clique_strat,
			ColorStrategy& color_strat,
			const std::string& log_file_path)
			: _branching_strat(branching_strat),
			_clique_strat(clique_strat),
			_color_strat(color_strat){
				_log_file.open(log_file_path);
				if (!_log_file.is_open()) {
					throw std::runtime_error("Failed to open log file: " + log_file_path);
				}
			}

		/**
         * @brief Solves the graph coloring problem using the branch and bound method.
         *
         * @param g The graph to solve.
         * @param optimum_time The time at which the optimum solution was found.
         * @param timeout_seconds The maximum time (in seconds) the solver can run before it times out.
         * @param sol_gather_period The period (in seconds) at which the best upper bound is gathered.
         * @return int The number of colors used in the optimal solution.
         */
		int Solve(Graph& g, double &optimum_time, int timeout_seconds = 60, int sol_gather_period = 10, unsigned short expected_chi = -1);
	};


class BalancedBranchNBoundPar {
	private:
		BranchingStrategy& _branching_strat;
		CliqueStrategy& _clique_strat;
		ColorStrategy& _color_strat;
		std::ofstream _log_file;

		std::atomic<unsigned short> _best_ub = USHRT_MAX;
		std::mutex _best_branch_mutex;
		Branch _current_best;

		void ColorInitialGraph(Graph& initial_graph, const Branch& optimal_branch);

		/**
		 * @brief updates the _current_best branch to store the best graph with the best coloring
		 * 
		 * @param depth depth at which the graph was colored
		 * @param lb lower bound of the best colored graph
		 * @param ub upper bound of the best colored graph
		 * @param graph the best colored graph
		 */
		void UpdateCurrentBest(int depth, int lb, unsigned short ub, GraphPtr graph);
	
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
		 * @param graph_to_color The graph to color.
		 */
		void thread_0_terminator(int my_rank, int p, int global_start_time, 
									int timeout_seconds, double &optimum_time,
									Graph& graph_to_color);
	
		/**
		 * @brief Updates (gathers) best_ub from time to time.
		 *
		 * @param p The total number of processes in the MPI communicator.
		 * @param sol_gather_period The period in seconds when all processes gather solutions.
		 */
		void thread_1_solution_gatherer(int p, int sol_gather_period);
	
		/**
		 * @brief Employer thread employs workers by answering their work requests.
		 *
		 * @param queue_mutex Mutex to protect concurrent access to the work queue.
		 * @param queue The local work queue containing branches to be processed.
		 */
		void thread_2_employer(std::mutex& queue_mutex, BranchQueue& queue);
	
	public:
		/**
		 * @brief Constructs a BalancedBranchNBoundPar solver.
		 *
		 * @param timeout_seconds The maximum time (in seconds) the solver can
		 * run before it times out.
		 * @param iteration_threshold The maximum number of iterations without
		 * improvement before stopping.
		 */
			BalancedBranchNBoundPar(BranchingStrategy& branching_strat,
			CliqueStrategy& clique_strat,
			ColorStrategy& color_strat,
			const std::string& log_file_path)
			: _branching_strat(branching_strat),
			_clique_strat(clique_strat),
			_color_strat(color_strat)
			{
				_log_file.open(log_file_path);
				if (!_log_file.is_open()) {
					throw std::runtime_error("Failed to open log file: " + log_file_path);
				}
			}
	
		int Solve(Graph& g, double &optimum_time, int timeout_seconds = 60, 
					int sol_gather_period = 10, 
					unsigned short expected_chi = -1);
	};
	

#endif	// BRANCH_N_BOUND_PAR_HPP
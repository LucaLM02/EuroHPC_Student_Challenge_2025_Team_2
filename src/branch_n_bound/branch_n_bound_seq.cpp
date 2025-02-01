#include "branch_n_bound_seq.hpp"

#include <algorithm>
#include <chrono>
#include <queue>
#include <utility>

bool BranchNBoundSeq::CheckTimeout(
    const std::chrono::steady_clock::time_point& start_time,
    int timeout_seconds) {
	auto current_time = std::chrono::steady_clock::now();
	auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
				   current_time - start_time)
				   .count();
	return elapsed_seconds >= timeout_seconds;
}

int BranchNBoundSeq::Solve(Graph& g, int timeout_seconds,
			   int iteration_threshold) {
	// Start the timeout timer
	auto start_time = std::chrono::steady_clock::now();

	// Initialize bounds
	int lb = _clique_strat->clique(g);
	int ub = _color_strat->color(g);
	int best_ub = ub;

	// Queue stores (current_G, current_lb, current_ub)
	std::queue<std::pair<Graph*, std::pair<int, int>>> queue;
	queue.push({&g, {lb, ub}});

	// Counter for iterations where best_ub remains unchanged
	int best_ub_unchanged_count = 0;

	while (!queue.empty()) {
		// Check for timeout
		if (CheckTimeout(start_time, timeout_seconds)) {
			std::cout << "Timeout reached. Returning best solution "
				     "found so far."
				  << std::endl;
			return best_ub;
		}
		// Dequeue the next branch
		auto current = queue.front();
		queue.pop();
		Graph* current_G = current.first;
		int current_lb = current.second.first;
		int current_ub = current.second.second;

		// If current_lb == current_ub, chromatic number found
		if (current_lb == current_ub) {
			return current_lb;
		}

		// Prune if current_lb >= best_ub
		if (current_lb >= best_ub) {
			continue;
		}

		// Find two non-adjacent vertices u and v according to strategy
		auto [u, v] = _branching_strat->ChooseVertices(*current_G);

		// If no such pair exists, the graph is complete (we are at a
		// leaf branch)
		if (u == -1 || v == -1) {
			return current_G->getNumVertices();
		}

		// Branch 1 - Merge u and v (assign same color)
		Graph* G1 = new Graph(*current_G);  // Copy Graph
		bool merged = G1->merge(u, v);
		if (merged) {
			int lb1 = _clique_strat->clique(*G1);
			int ub1 = _color_strat->color(*G1);
			if (lb1 < best_ub) {
				queue.push({G1, {lb1, std::min(ub1, best_ub)}});
			}
		}

		// Branch 2 - Add edge between u and v (assign different colors)
		Graph* G2 = new Graph(*current_G);  // Copy Graph
		G2->addEdge(u, v);
		int lb2 = _clique_strat->clique(*G2);
		int ub2 = _color_strat->color(*G2);
		if (lb2 < best_ub) {
			queue.push({G2, {lb2, std::min(ub2, best_ub)}});
		}

		// Update best_ub
		int previous_best_ub = best_ub;
		best_ub = std::min({best_ub, ub1, ub2});

		// Check if best_ub has changed
		if (best_ub == previous_best_ub) {
			best_ub_unchanged_count++;
		} else {
			best_ub_unchanged_count = 0;  // Reset the counter
		}

		// Check if the iteration threshold is exceeded
		if (best_ub_unchanged_count >= iteration_threshold) {
			std::cout << "Iteration threshold reached. Returning "
				     "best solution found so far."
				  << std::endl;
			return best_ub;
		}
	}

	return best_ub;
}

#include "branch_n_bound_seq.hpp"

bool BranchNBoundSeq::CheckTimeout(
    const std::chrono::steady_clock::time_point& start_time,
    int timeout_seconds) {
	auto current_time = std::chrono::steady_clock::now();
	auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
				   current_time - start_time)
				   .count();
	return elapsed_seconds >= timeout_seconds;
}

void BranchNBoundSeq::Log(const std::string& message) {
	if (_log_file.is_open()) {
		_log_file << message << std::endl;
	}
}

int BranchNBoundSeq::Solve(Graph& g, int timeout_seconds,
			   int iteration_threshold) {
	// Start the timeout timer
	auto start_time = std::chrono::steady_clock::now();

	// Initialize bounds
	int lb = _clique_strat->Clique(g);
	int ub = _color_strat->Color(g);
	int best_ub = ub;

	// Log initial bounds
	Log("Initial bounds: lb = " + std::to_string(lb) +
	    ", ub = " + std::to_string(ub));

	// Queue stores (current_G, current_lb, current_ub)
	std::queue<std::pair<Graph*, std::pair<int, int>>> queue;
	queue.push({&g, {lb, ub}});

	// Counter for iterations where best_ub remains unchanged
	int best_ub_unchanged_count = 0;

	while (!queue.empty()) {
		// Check for timeout
		if (CheckTimeout(start_time, timeout_seconds)) {
			Log("Timeout reached. Returning best solution found so "
			    "far: " +
			    std::to_string(best_ub));
			return best_ub;
		}

		// Dequeue the next branch
		auto current = queue.front();
		queue.pop();
		Graph* current_G = current.first;
		int current_lb = current.second.first;
		int current_ub = current.second.second;

		// Log current node
		Log("Processing node: lb = " + std::to_string(current_lb) +
		    ", ub = " + std::to_string(current_ub));

		// If current_lb == current_ub, chromatic number found
		if (current_lb == current_ub) {
			Log("Chromatic number found: " +
			    std::to_string(current_lb));
			return current_lb;
		}

		// Prune if current_lb >= best_ub
		if (current_lb >= best_ub) {
			Log("Pruning branch: lb = " +
			    std::to_string(current_lb) +
			    " >= best_ub = " + std::to_string(best_ub));
			continue;
		}

		// Find two non-adjacent vertices u and v according to strategy
		auto [u, v] = _branching_strat->ChooseVertices(*current_G);
		Log("Branching on vertices: u = " + std::to_string(u) +
		    ", v = " + std::to_string(v));

		// If no such pair exists, the graph is complete (we are at a
		// leaf branch)
		if (u == -1 || v == -1) {
			Log("Graph is complete. Chromatic number = " +
			    std::to_string(current_G->getNumVertices()));
			return current_G->getNumVertices();
		}

		// Branch 1 - Merge u and v (assign same color)
		Graph* G1 = new Graph(*current_G);  // Copy Graph
		bool merged = G1->merge(u, v);
		if (merged) {
			int lb1 = _clique_strat->clique(*G1);
			int ub1 = _color_strat->color(*G1);
			Log("Branch 1 (merge): lb = " + std::to_string(lb1) +
			    ", ub = " + std::to_string(ub1));
			if (lb1 < best_ub) {
				queue.push({G1, {lb1, ub1}});
			}
		}

		// Branch 2 - Add edge between u and v (assign different colors)
		Graph* G2 = new Graph(*current_G);  // Copy Graph
		G2->addEdge(u, v);
		int lb2 = _clique_strat->clique(*G2);
		int ub2 = _color_strat->color(*G2);
		Log("Branch 2 (add edge): lb = " + std::to_string(lb2) +
		    ", ub = " + std::to_string(ub2));
		if ((lb2 < best_ub) && (lb2 < ub1)) {
			queue.push({G2, {lb2, ub2}});
		}

		// Update best_ub
		int previous_best_ub = best_ub;
		best_ub = std::min({best_ub, ub1, ub2});
		Log("Updated best_ub: " + std::to_string(best_ub));

		// Check if best_ub has changed
		if (best_ub == previous_best_ub) {
			best_ub_unchanged_count++;
		} else {
			best_ub_unchanged_count = 0;  // Reset the counter
		}

		// Check if the iteration threshold is exceeded
		if (best_ub_unchanged_count >= iteration_threshold) {
			Log("Iteration threshold reached. Returning best "
			    "solution found so far: " +
			    std::to_string(best_ub));
			return best_ub;
		}
	}

	Log("Final chromatic number: " + std::to_string(best_ub));
	return best_ub;
}

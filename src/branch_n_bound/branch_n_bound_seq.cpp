#include "branch_n_bound_seq.hpp"

bool BranchNBoundPar::CheckTimeout(
    const std::chrono::steady_clock::time_point& start_time,
    int timeout_seconds) {
	auto current_time = MPI_Wtime();
	auto elapsed_seconds = current_time - start_time;
	return elapsed_seconds >= timeout_seconds;
}

void BranchNBoundPar::Log(const std::string& message) {
	if (_log_file.is_open()) {
		_log_file << message << std::endl;
	}
}

void terminate_graph_proccessing() {
	int termination_flag[1];
	termination_flag[0] = 1;
	MPI_Bcast(termination_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// TODO: Should the input be a vector of graphs or should we call Solve in real
// time continuously by popping from a queue?
int BranchNBoundPar::Solve(Graph& g, int timeout_seconds,
			   int iteration_threshold) {
	// Intitialize MPI
	MPI_Init(NULL, NULL);
	int my_rank;
	int p;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Start the timeout timer
	auto start_time = MPI_Wtime();

	// TODO: Implement outer loop to process multiple graphs
	if (my_rank == 0) {
		// Master initialize bounds
		int lb = _clique_strat->clique(g);
		int ub = _color_strat->color(g);
		int best_ub = ub;

		// Log initial bounds
		Log("Initial bounds: lb = " + std::to_string(lb) +
		    ", ub = " + std::to_string(ub));

		// Priority queue stores (depth, (current_G, current_lb,
		// current_ub))
		std::priority_queue<struct Branch> p_queue;
		struct Branch root = {g*, lb, ub, 0};
		p_queue.push(root);

		// Counter for iterations where best_ub remains unchanged
		int best_ub_unchanged_count = 0;

		while (!queue.empty()) {
			// Check for timeout
			if (CheckTimeout(start_time, timeout_seconds)) {
				Log("Timeout reached. Returning best solution "
				    "found so "
				    "far: " +
				    std::to_string(best_ub));
				return best_ub;
			}

			// Dequeue the deepest branch
			auto current = queue.pop();
			Graph* current_G = current.g;
			int current_lb = current.lb;
			int current_ub = current.ub;

			// Log current node
			Log("Processing node: lb = " +
			    std::to_string(current_lb) +
			    ", ub = " + std::to_string(current_ub));

			// If current_lb == current_ub, chromatic number found
			if (current_lb == current_ub) {
				Log("Chromatic number found: " +
				    std::to_string(current_lb));
				best_ub = current_lb;
				terminate_graph_proccessing();
				break;
			}

			// Prune if current_lb >= best_ub
			if (current_lb >= best_ub) {
				Log("Pruning branch: lb = " +
				    std::to_string(current_lb) +
				    " >= best_ub = " + std::to_string(best_ub));
				continue;
			}

			// Find two non-adjacent vertices u and v according to
			// strategy
			auto [u, v] =
			    _branching_strat->ChooseVertices(*current_G);
			Log("Branching on vertices: u = " + std::to_string(u) +
			    ", v = " + std::to_string(v));

			// If no such pair exists, the graph is complete (we are
			// at a leaf branch)
			if (u == -1 || v == -1) {
				Log("Graph is complete. Chromatic number = " +
				    std::to_string(
					current_G->getNumVertices()));
				best_ub = current_G->getNumVertices();
				terminate_graph_processing();
				break;
			}

			// Branch 1 - Merge u and v (assign same color)
			Graph* G1(current_G);  // Copy Graph
			bool merged = G1->merge(u, v);
			if (merged) {
				int available_rank;
				MPI_Status stat;
				int MPI_Recv(&available_rank, 1, MPI_INT,
					     MPI_ANY_SOURCE, MPI_ANY_TAG,
					     MPI_COMM_WORLD, &stat);

				// TODO: MPI Graph subclasses Derived Types
				// should be defined.
				// TODO: Here send Graph object lengths and
				// object itself (using the derived type).

				// TODO: clique and color should be moved
				// outside master if
				// int lb1 = _clique_strat->clique(*G1);
				// int ub1 = _color_strat->color(*G1);
				Log("Branch 1 (merge): lb = " +
				    std::to_string(lb1) +
				    ", ub = " + std::to_string(ub1));
				if (lb1 < best_ub) {
					queue.push({G1, {lb1, ub1}});
				}
			}

			// Branch 2 - Add edge between u and v (assign different
			// colors)
			Graph* G2(current_G);  // Copy Graph
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
				best_ub_unchanged_count =
				    0;	// Reset the counter
			}

			// Check if the iteration threshold is exceeded
			if (best_ub_unchanged_count >= iteration_threshold) {
				Log("Iteration threshold reached. Returning "
				    "best "
				    "solution found so far: " +
				    std::to_string(best_ub));
				return best_ub;
			}
		}

		if (my_rank == 0)
			Log("Final chromatic number: " +
			    std::to_string(best_ub));
		MPI_Finalize();
		return best_ub;
	}

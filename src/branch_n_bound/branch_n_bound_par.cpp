#include "branch_n_bound_par.hpp"

#define ALLGATHER_WAIT_TIME 10	   // Sleep time for MPI_Allgather
#define TIMEOUT_CHECK_WAIT_TIME 5  // Sleep time for timeout checker

//tag for MPI communication
#define TAG_WORK_REQUEST 1
#define TAG_WORK_RESPONSE 2
#define TAG_WORK 3
#define TAG_SOLUTION_FOUND 4

//now i fixed with start time seconds, check if its better using only mpi_wtime
bool BranchNBoundPar::CheckTimeout( 
    const std::chrono::steady_clock::time_point& start_time,
    int timeout_seconds) {
	auto current_time = MPI_Wtime();
	double start_time_seconds = std::chrono::duration<double>(start_time.time_since_epoch()).count();
	auto elapsed_seconds = current_time - start_time_seconds;
	return elapsed_seconds >= timeout_seconds;
}

void BranchNBoundPar::Log(const std::string& message, int depth = 0,
			  bool is_branching = false) {
	if (_log_file.is_open()) {
		// Indentation based on depth
		std::string indentation(depth * 2, ' ');

		if (is_branching) {
			_log_file << indentation << "[Depth " << depth
				  << "] Branching on u = " << message
				  << std::endl;
		} else {
			_log_file << indentation << message << std::endl;
		}
	}
}

/**
 * thread_0_solution_gatherer - Periodically gathers the best upper bound
 * (best_ub) from all worker processes and updates the global best_ub. This
 * function is called by the first thread of the master and the worker
 * processes.
 *
 * This function is run by the master process and performs an allgather
 * operation to collect the best upper bound from all worker nodes every
 * ALLGATHER_WAIT_TIME seconds. The best upper bound is updated by taking the
 * minimum of the current best_ub and the gathered values.
 *
 * Parameters:
 *   p (int)           : The number of processes in the MPI communicator.
 *   best_ub (int*)    : Pointer to the variable holding the best upper bound.
 */
void thread_0_solution_gatherer(int p, unsigned short& best_ub) {
	std::vector<int> all_best_ub(p);
	// Start timer for best_ub exchange interval.
	time_t start_time = time(NULL);
	while (1) {
		// Wait for the time threshold. Allgather in done every
		// ALLGATHER_WAIT_TIME seconds.
		sleep(ALLGATHER_WAIT_TIME);
		// Gather best_ub from other workers
		MPI_Allgather(&best_ub, 1, MPI_INT, all_best_ub.data(), 1, MPI_INT,
			      MPI_COMM_WORLD);
		// Update the best upper bound for other threads in
		// shared memory
		best_ub = *std::min_element(all_best_ub.begin(), all_best_ub.end());
		// Reset timer
		start_time = time(NULL);
	}
}

/**
 * thread_1_terminator - Listens for termination signals (solution found or
 * timeout) and terminates the execution of the process accordingly. This
 * function is called by the second thread of the master and the worker
 * processes.
 *
 * This function is responsible for detecting whether a solution has been found
 * or a timeout has occurred. If the master (rank 0) receives a solution, it
 * broadcasts the solution status to all processes. If the timeout occurs, it
 * broadcasts the timeout signal. Worker nodes listen for these termination
 * signals and exit the loop if any of the conditions are met.
 *
 * Parameters:
 *   my_rank (int)           : The rank of the current process.
 *   p (int)                 : The total number of processes in the MPI
 * communicator. global_start_time (int) : The global start time, used to check
 * for timeouts. timeout_seconds (int)   : The timeout duration (in seconds)
 * after which the timeout signal is sent.
 */
void thread_1_terminator(int my_rank, int p, int global_start_time,
			 int timeout_seconds) {
	int solution_found = 0;
	int timeout_signal = 0;
	while (1) {
		if (my_rank == 0) {
			// Master listens for solution found (Non-blocking)
			MPI_Status status;
			int flag;
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
				   &flag, &status);
			// Check if a solution is being communicated
			if (flag) {
				int solution;
				MPI_Recv(&solution, 1, MPI_INT,
					 status.MPI_SOURCE, status.MPI_TAG,
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				solution_found = 1;
				// Notify all processes that a solution has been
				// found
				MPI_Bcast(&solution_found, 1, MPI_INT, 0,
					  MPI_COMM_WORLD);
				break;	// Exit loop after broadcasting solution
					// found
			}

			// Check if timeout is reached, broadcast timeout signal
			// TODO: Do we really care if workers exit forcefully
			// with MPI finalize? If we dont, we dont need all this
			// broadcasting stuff.
			if (MPI_Wtime() - global_start_time >=
			    timeout_seconds) {
				timeout_signal = 1;
				MPI_Bcast(&timeout_signal, 1, MPI_INT, 0,
					  MPI_COMM_WORLD);  // Broadcast timeout
							    // signal
				break;	// Exit loop after broadcasting timeout
					// signal
			}

		} else {
			// Worker nodes listen for termination signals (solution
			// or timeout)
			MPI_Bcast(&solution_found, 1, MPI_INT, 0,
				  MPI_COMM_WORLD);  // Receive solution status
						    // from master
			MPI_Bcast(&timeout_signal, 1, MPI_INT, 0,
				  MPI_COMM_WORLD);  // Receive timeout status
						    // from master

			// Exit if solution or timeout is detected
			if (solution_found || timeout_signal) {
				break;	// Exit the loop if solution or timeout
					// has been detected
			}
		}
		usleep(10000);	// Prevent CPU overload (10 ms)
	}
}

// Function to listen for requests from other workers.
//TODO: Implement this function as response if has work than send work (reference line 388)
void thread_2_listen_for_requests() {
	int request_signal = 0;
	while (true) {
		// Listen for a request for work from other workers.
		MPI_Bcast(&request_signal, 1, MPI_INT, 0,
			  MPI_COMM_WORLD);  // Listen for request signal

		if (request_signal == 1 && !queue.empty()) {
			// If a request is received and this worker has work,
			// send work.

			//need mutex to avoid cuncurrent access
			Branch branch = std::move(queue.top());
			queue.pop();
			// Send branch to the requesting worker (using MPI
			// Send).
			int destination_rank =
			    1;	// This should be the rank of the requesting
				// worker, for simplicity assume rank 1
			MPI_Send(&branch, sizeof(Branch), MPI_BYTE,
				 destination_rank, 0, MPI_COMM_WORLD);
		}

		usleep(10000);	// Prevent CPU overload (10 ms)
	}
}

int BranchNBoundPar::Solve(Graph& g, int timeout_seconds, int iteration_threshold) {
	// Intitialize MPI
	MPI_Init(NULL, NULL);
	int my_rank;
	int p;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Start the timeout timer
	auto global_start_time = MPI_Wtime();
	// Initialize big enough best_ub for all processes.
	unsigned short best_ub = USHRT_MAX;

	std::mutex queue_mutex; //avoid concurrent access to the queue
	using BranchQueue = std::priority_queue<Branch, std::vector<Branch>>;
	BranchQueue queue;

	if (my_rank == 0) {
		// Initialize bounds
		int lb = _clique_strat.FindClique(g);
		unsigned short ub;
		_color_strat.Color(g, ub);
		best_ub = ub;

		// Log initial bounds
		Log("[INIT] Initial bounds: lb = " + std::to_string(lb) +
		    ", ub = " + std::to_string(ub));

		queue.push(Branch(g.Clone(), lb, ub,
				  1));	// Initial branch with depth 1

		// Start computation tree
		Log("========== START ==========", 0);
		Log("Root: (lb = " + std::to_string(lb) +
			", ub = " + std::to_string(ub) + ")",
		    0);

		// Get enough branches to distribute amongst
		// workers
		while (queue.size() < p - 1) {
			// Dequeue the next branch
			Branch current =
			    std::move(const_cast<Branch&>(queue.top()));
			queue.pop();
			auto current_G = std::move(current.g);
			int current_lb = current.lb;
			int current_ub = current.ub;

			// Log current node
			Log("Processing node: lb = " +
				std::to_string(current_lb) +
				", ub = " + std::to_string(current_ub),
			    current.depth);

			// If current_lb == current_ub,
			// chromatic number found
			if (current_lb == current_ub) {
				Log("[FOUND] Chromatic number "
				    "found: " +
					std::to_string(current_lb),
				    current.depth);
				Log("========== END ==========", 0);
				best_ub = current_lb;
				break;
			}

			// Prune if current_lb >= best_ub
			if (current_lb >= best_ub) {
				Log("[PRUNE] Branch pruned at "
				    "depth " +
					std::to_string(current.depth) +
					": lb = " + std::to_string(current_lb) +
					" >= best_ub = " +
					std::to_string(best_ub),
				    current.depth);
				continue;
			}

			// Find two non-adjacent vertices u and
			// v according to strategy
			auto type = _branching_strat.PairType::DontCare;
			auto [u, v] =
			    _branching_strat.ChooseVertices(*current_G, type);
			Log("Branching on vertices: u = " + std::to_string(u) +
				", v = " + std::to_string(v),
			    current.depth, true);

			// If no such pair exists, the graph is
			// complete (we are at a leaf branch)
			if (u == -1 || v == -1) {
				best_ub = current_G->GetNumVertices();
				Log("Graph is complete. "
				    "Chromatic number = " +
					std::to_string(best_ub),
				    current.depth);
				Log("========== END ==========", 0);
				break;
			}

			// Branch 1 - Merge u and v (assign same
			// color)
			auto G1 = current_G->Clone();  // Copy Graph
			G1->MergeVertices(u, v);
			int lb1 = _clique_strat.FindClique(*G1);
			unsigned short ub1;
			_color_strat.Color(*G1, ub1);
			Log("[Branch 1] (Merge u, v) lb = " +
				std::to_string(lb1) +
				", ub = " + std::to_string(ub1),
			    current.depth);
			if (lb1 < best_ub) {
				queue.push(Branch(std::move(G1), lb1, ub1,
						  current.depth + 1));
			}

			// Branch 2 - Add edge between u and v
			// (assign different colors)
			auto G2 = current_G->Clone();  // Copy Graph
			G2->AddEdge(u, v);
			int lb2 = _clique_strat.FindClique(*G2);
			unsigned short ub2;
			_color_strat.Color(*G2, ub2);
			Log("[Branch 2] (Add edge u-v) lb = " +
				std::to_string(lb2) +
				", ub = " + std::to_string(ub2),
			    current.depth);
			if ((lb2 < best_ub) && (lb2 < ub1)) {
				queue.push(Branch(std::move(G2), lb2, ub2,
						  current.depth + 1));
			}

			// Update best_ub
			best_ub = std::min({best_ub, ub1, ub2});
			Log("[UPDATE] Updated best_ub: " +
				std::to_string(best_ub),
			    current.depth);
		}
		// Send all p branches to all p workers.
		while (!queue.empty()) {
			Branch current =
			    std::move(const_cast<Branch&>(queue.top()));
			queue.pop();
			auto current_G = std::move(current.g);
			int current_lb = current.lb;
			int current_ub = current.ub;
			// TODO: Find an efficient and abstract way to send
			// different Graph types to other workers.
			// TODO: Send popped graph to a worker.
		}
		Log("[PARALLELISATION START]");
	}

	if (my_rank != 0) {
		// TODO: Here create priority queue for branches.
		// TODO: Here receive graph from Master, compute bounds and push
		// it to priority queue (G, lb, ub).

		// TODO: implement log

		MPI_Status status_recv;
		Branch branch_recv; //initialize branch
    
    	// request work
    	MPI_Send(nullptr, 0, MPI_INT, 0, TAG_WORK_REQUEST, MPI_COMM_WORLD);
    
    	// recv work
    	MPI_Recv(&branch_recv, sizeof(Branch), MPI_BYTE, 0, TAG_WORK, MPI_COMM_WORLD, &status_recv); //send branch
		unsigned short best_ub = branch_recv.ub;

		queue.push(std::move(branch_recv));
		
		int best_ub_unchanged_count = 0; //do we need this?
	}

	int active_tasks = 0;
	int max_tasks = omp_get_max_threads()-3; //limit number of tasks (maybe we can increase this number)

	// OpenMP Parallel Region
	#pragma omp parallel shared(best_ub, queue, queue_mutex, active_tasks, max_tasks){
			MPI_Status status_work_steal;
			int tid = omp_get_thread_num();
			// Both master's and worker's thread 0 goes in here.
			if (tid == 0) {
				// Updates (gathers) best_ub from time to time.
				thread_0_solution_gatherer(p, best_ub);
			}
			// Both master's and worker's thread 1 goes in here.
			if (tid == 1) {
				// Checks if solution has been found or timeout.
				thread_1_terminator(my_rank, p, global_start_time, timeout_seconds);
			}
			//worker thread used for listening for requests about work stealing
			if (tid == 2 && my_rank != 0) {
				thread_2_listen_for_requests();
			}
			// Only worker processes go in here.
			if (my_rank != 0) {
				// TODO: Here worker should pop first element and start
				// branching as well as implement work stealing (maybe
				// create own communication group for this?)

				Branch current;
				std::unique_ptr<Graph> current_G;
				int current_lb;
				unsigned short current_ub;
				unsigned int u,v;
				std::vector<Branch> new_branches(2);

				#pragma omp single nowait shared(queue, queue_mutex, best_ub, active_tasks, max_tasks){

					//TODO:change true with flag terminate (is false until we receive a signal from terminator)
					while (true) {  // keep dequeuing until queue is empty or we dont have other threads available
						bool has_work = false;

						if (!queue.empty()) {
							std::lock_guard<std::mutex> lock(queue_mutex);
							current = std::move(const_cast<Branch&>(queue.top()));
							queue.pop();
							has_work = true;
						}

						//TODO: if response is 0 ask to others(limit number of requests?)
						if (!has_work) { // Work Stealing

							int target_worker = (my_rank + 1) % p; //check this

							//ask for work
							int response = 0;
							MPI_Send(nullptr, 0, MPI_INT, target_worker, TAG_WORK_REQUEST, MPI_COMM_WORLD);
							MPI_Recv(&response, 1, MPI_INT, target_worker, TAG_WORK_RESPONSE, MPI_COMM_WORLD, &status_work_steal);

							if (response == 1) { //there is work
								MPI_Recv(&current, sizeof(Branch), MPI_BYTE, target_worker, TAG_WORK, MPI_COMM_WORLD, &status_work_steal);
								std::lock_guard<std::mutex> lock(queue_mutex);
								queue.push(std::move(current));
							} else {
								//no work available
								continue;
							}
						}

						current_G = std::move(current.g);
						current_lb = current.lb;
						current_ub = current.ub;

						// find solution
						if (current_lb == current_ub) {
							MPI_Send(&current_ub, 1, MPI_INT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD); //we need also coloring?
							continue;
						}

						// Prune
						if (current_lb >= best_ub) continue;

						auto type = _branching_strat.PairType::DontCare;
						std::pair<unsigned int, unsigned int> vertices = _branching_strat.ChooseVertices(*current_G, type);
						u = vertices.first;
						v = vertices.second;

						//branching
						if (u == -1 || v == -1) {
							int chromatic_number = current_G->GetNumVertices();
							MPI_Send(&chromatic_number, 1, MPI_INT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD); //we need also coloring?
							continue;
						}

						while (active_tasks >= max_tasks) {
							#pragma omp taskyield // avoid creating too many tasks
						}

						#pragma omp critical //use atomic?
						active_tasks = active_tasks + 2;

						#pragma omp taskgroup shared(new_branches, current_G, u, v)
						{

							#pragma omp task firstprivate(current_G, u, v)
							{
								auto G1 = current_G->Clone();
								G1->MergeVertices(u, v);
								int lb1 = _clique_strat.FindClique(*G1);
								unsigned short ub1;
								_color_strat.Color(*G1, ub1);

								if (lb1 < best_ub) {
									new_branches[0] = Branch(std::move(G1), lb1, ub1, current.depth + 1);
								}
							}

							#pragma omp task firstprivate(current_G, u, v)
							{
								auto G2 = current_G->Clone();
								G2->AddEdge(u, v);
								int lb2 = _clique_strat.FindClique(*G2);
								unsigned short ub2;
								_color_strat.Color(*G2, ub2);

								if ((lb2 < best_ub) && (lb2 < new_branches[0].ub)) { //possible problem with condition (we need to check if ub1 is defined)
									new_branches[1] = Branch(std::move(G2), lb2, ub2, current.depth + 1);
								}
							}
						
						}

						#pragma omp critical //use atomic?
						active_tasks= active_tasks - 2;

						//add new branches to the queue
						{
							std::lock_guard<std::mutex> lock(queue_mutex);
							for (auto& branch : new_branches) {
								if (branch.g) queue.push(std::move(branch));
							}
					    }

						int previous_best_ub = best_ub;
						#pragma omp critical
						best_ub = std::min({best_ub, new_branches[0].ub, new_branches[1].ub});
					}
				}
			}
		}
		// End execution
		if (my_rank == 0){
			Log("Final chromatic number: " + std::to_string(best_ub));
		}
		MPI_Finalize();
		return best_ub;
}

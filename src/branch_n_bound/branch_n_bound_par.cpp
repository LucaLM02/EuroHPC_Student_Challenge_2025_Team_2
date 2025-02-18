#include "branch_n_bound_par.hpp"

#define ALLGATHER_WAIT_TIME 5	   // Sleep time for MPI_Allgather
#define TIMEOUT_CHECK_WAIT_TIME 5  // Sleep time for timeout checker

// tags for MPI communication
#define TAG_WORK_REQUEST 1
#define TAG_WORK_RESPONSE 2
#define TAG_WORK 3
#define TAG_SOLUTION_FOUND 4

using BranchQueue = std::priority_queue<Branch, std::vector<Branch>>;
std::atomic<bool> terminate_flag = false;
std::atomic<int> active_tasks = 0;
int max_tasks = 10;	// limit number of tasks (maybe we can increase this number)

/**
 * @brief Checks if the timeout has been reached based on the start time and the
 * given timeout duration.
 *
 * @param start_time The start time point from which to calculate the elapsed
 * time.
 * @param timeout_seconds The timeout duration in seconds.
 * @return true if the elapsed time is greater than or equal to the timeout,
 * false otherwise.
 */
bool BranchNBoundPar::CheckTimeout(
    const std::chrono::steady_clock::time_point& start_time,
    int timeout_seconds) {
	// TODO: check if its better using only mpi_wtime instead of chrono
	auto current_time = MPI_Wtime();
	double start_time_seconds =
	    std::chrono::duration<double>(start_time.time_since_epoch())
		.count();
	auto elapsed_seconds = current_time - start_time_seconds;
	return elapsed_seconds >= timeout_seconds;
}

void BranchNBoundPar::Log(const std::string& message, int depth = 0,
			  bool is_branching = false) {
	if (_log_file_master.is_open()) {
		// Indentation based on depth
		std::string indentation(depth * 2, ' ');

		if (is_branching) {
			_log_file_master << indentation << "[Depth " << depth
				  << "] Branching on u = " << message
				  << std::endl;
		} else {
			_log_file_master << indentation << message << std::endl;
		}
	}
}

void BranchNBoundPar::Log_par(const std::string& message, int depth = 0, bool is_branching = false) {
	if (_log_file_branches.is_open()) {
		int rank = 0, thread_id = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		thread_id = omp_get_thread_num();

		// Indentation based on depth
		std::string indentation(depth * 2, ' ');
		// TODO: Logging everything with this critical section is not efficient
		#pragma omp critical
		{
			if (is_branching) {
				_log_file_branches << indentation << "[Rank " << rank
					  << " | Thread " << thread_id << "] "
					  << "[Depth " << depth
					  << "] Branching on u = " << message
					  << std::endl;
			} else {
				_log_file_branches << indentation << "[Rank " << rank
					  << " | Thread " << thread_id << "] "
					  << message << std::endl;
			}
		}
	}
}

/**
 * @brief Sends a serialized Branch object to a specified destination in the MPI
 * communication.
 *
 * This function serializes the Branch object and sends it over MPI to the
 * specified destination with the provided tag and communicator.
 *
 * @param b The Branch object to be sent.
 * @param dest The destination rank for the message.
 * @param tag The MPI tag for the message.
 * @param comm The MPI communicator used for communication.
 */
void sendBranch(const Branch& b, int dest, int tag, MPI_Comm comm) {
	//std::cout << "Before serializing branch" << std::endl;
	std::vector<char> buffer = b.serialize();
	int size = buffer.size();
	//std::cout << "Sucessfully serialized branch with size " << size << " bytes."
	//		  << std::endl;

	MPI_Send(&size, 1, MPI_INT, dest, tag, comm);
	MPI_Send(buffer.data(), size, MPI_BYTE, dest, tag, comm);
}

/**
 * @brief Receives a serialized Branch object from a specified source in the MPI
 * communication.
 *
 * This function receives the serialized Branch object from the specified source
 * and reconstructs the Branch object by deserializing the received data.
 *
 * @param source The source rank from which to receive the message.
 * @param tag The MPI tag for the message.
 * @param comm The MPI communicator used for communication.
 * @return The deserialized Branch object.
 */
Branch recvBranch(int source, int tag, MPI_Comm comm) {
	MPI_Status status;
	int size;

	MPI_Recv(&size, 1, MPI_INT, source, tag, comm, &status);

	std::vector<char> buffer(size);
	MPI_Recv(buffer.data(), size, MPI_BYTE, source, tag, comm, &status);
	//std::cout << "Sucessfully received serialized branch with size " << size << " bytes."
	//	  		<< std::endl;
	return Branch::deserialize(buffer);
}

/**
 * thread_0_terminator - Listens for termination signals (solution found or
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
void thread_0_terminator(int my_rank, int p, int global_start_time, int timeout_seconds) {
	int solution_found = 0;
	int timeout_signal = 0;

	while (true) {
		if (my_rank == 0) {
			// Master listens for solution found (Non-blocking)
			MPI_Status status;
			int flag = 0;
			MPI_Iprobe(MPI_ANY_SOURCE, TAG_SOLUTION_FOUND, MPI_COMM_WORLD, &flag, &status);
			// Check if a solution is being communicated
			if (flag) {
				unsigned short solution;
				MPI_Recv(&solution, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				solution_found = 1;
				std::cout << "Solution found by worker: " << status.MPI_SOURCE << std::endl;
				std::cout << "Solution: " << solution << std::endl;
			}
			// Check if timeout is reached, broadcast timeout signal
			// TODO: Do we really care if workers exit forcefully with MPI finalize? If we dont, we dont need all this broadcasting stuff
			if (MPI_Wtime() - global_start_time >= timeout_seconds) {
				std::cout << "Rank: " << " Timeout!" << std::endl;
				timeout_signal = 1;
			}
		}
		// Worker nodes listen for termination signals (solution or timeout)
		MPI_Bcast(&solution_found, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&timeout_signal, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Exit if solution or timeout is detected
		if (solution_found || timeout_signal) {
			terminate_flag.store(true, std::memory_order_relaxed);
			break;
		}
		usleep(10000);	// Prevent CPU overload (10 ms)
	}
}

/**
 * thread_1_solution_gatherer - Periodically gathers the best upper bound
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

// CHANGE: all_best_ub is now a vector of unsigned short, best_ub is now an
// atomic variable(helps to avoid cuncurrent access)
void thread_1_solution_gatherer(int p, std::atomic<unsigned short>& best_ub, int& my_rank) {
    std::vector<unsigned short> all_best_ub(p);
    auto last_gather_time = std::chrono::steady_clock::now();
    while (!terminate_flag.load(std::memory_order_relaxed)) {
		//std::cout << "prova2" << std::endl;
        auto current_time = std::chrono::steady_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(current_time - last_gather_time).count();
        if (elapsed_time >= ALLGATHER_WAIT_TIME) {
            unsigned short local_best_ub = best_ub.load();	// safe read
            // Gather best_ub from other workers
            MPI_Allgather(&local_best_ub, 1, MPI_UNSIGNED_SHORT, all_best_ub.data(), 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
            // Update the best upper bound for other threads in shared memory
            best_ub.store(*std::min_element(all_best_ub.begin(), all_best_ub.end()));  // safe write
            // Reset the timer
            last_gather_time = current_time;
        }
        // Sleep for a short duration to prevent busy-waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

// Function to listen for requests from other workers.
// TODO: evaluate using iprobe instead of recv to reduce comunication overhead
// CHANGE: added mutex to avoid concurrent access to the queue, added response
// to avoid sending work to workers if there is no work available(avoid also
// deadlock in workers)
void thread_2_listen_for_requests(std::mutex& queue_mutex, BranchQueue& queue,
				  int& my_rank) {
	MPI_Status status;
	int request_signal;
	while (!terminate_flag.load(std::memory_order_relaxed)) {

		// Listen for a request for work from other workers.
		MPI_Iprobe(MPI_ANY_SOURCE, TAG_WORK_REQUEST, MPI_COMM_WORLD, &request_signal, &status);
		//std::cout << "Rank: " << my_rank << " Looping in thread_2_listen_for_requests" << std::endl;
		int destination_rank = status.MPI_SOURCE;
		int response = 0;
		// Check if a solution is being communicated
		if (request_signal) {
			std::cout << "Rank: " << my_rank << " Received work request from rank " << destination_rank << std::endl;
			{
				std::lock_guard<std::mutex> lock(queue_mutex);
				if (!queue.empty()) {
					response = 1;
					MPI_Send(&response, 1, MPI_INT, destination_rank, TAG_WORK_REQUEST, MPI_COMM_WORLD);
						std::cout << "Rank: " << my_rank << " Responding positively to rank " << destination_rank << std::endl;
					Branch branch = std::move(const_cast<Branch&>(queue.top()));
					queue.pop();

					MPI_Send(&branch, sizeof(Branch), MPI_BYTE, destination_rank, TAG_WORK, MPI_COMM_WORLD);
					std::cout << "Rank: " << my_rank << " Sent work to rank " << destination_rank << std::endl;
				} else {
					std::cout << "Rank: " << my_rank << " Responding negatively to rank " << destination_rank << std::endl;
					MPI_Send(&response, 1, MPI_INT, destination_rank,
						TAG_WORK_REQUEST, MPI_COMM_WORLD);
				}
			}
		}	
	}
}

int work_stealing(int my_rank, int p, BranchQueue& queue, std::mutex& queue_mutex, Branch& current) {
	int target_worker = (my_rank + rand() % (p - 1) + 1) % p;	// check if it is correct
	MPI_Status status;
	// ask for work
	int response = 0;
	MPI_Send(nullptr, 0, MPI_INT, target_worker, TAG_WORK_REQUEST, MPI_COMM_WORLD);
	
	MPI_Iprobe(target_worker, TAG_WORK_RESPONSE, MPI_COMM_WORLD, &response, &status);
	if (response == 1) {  // there is work
		std::cout << "Rank: " << my_rank << " positive work stealing response from " << target_worker <<std::endl;
		current = recvBranch(status.MPI_SOURCE, TAG_WORK, MPI_COMM_WORLD);
		std::cout << "Rank: " << my_rank << " received work from " << target_worker <<std::endl;
		std::lock_guard<std::mutex> lock(queue_mutex);
		queue.push(std::move(current));	
	}
	return response;
}

/**
 * Creates a task for the OpenMP parallel region to execute.
 *
 *
 * @param u The first vertex to merge or add an edge.
 * @param v The second vertex to merge or add an edge.
 * @param _clique_strat The clique strategy.
 * @param _color_strat The color strategy.
 * @param new_branches The new branches to add to the queue.
 * @param task_type The type of task to create (1 = MergeVertices, 2 = AddEdge).
 * @param best_ub The best upper bound.
 */
void BranchNBoundPar::create_task(
    std::unique_ptr<Graph> current_G, int u, int v,
    CliqueStrategy& _clique_strat, ColorStrategy& _color_strat,
	std::atomic<unsigned short> & best_ub, int const & depth, int my_rank) {
/*

		//std::cout << "Rank: " << my_rank << " start task" << std::endl;
		active_tasks.fetch_add(1);
		std::cout << "Rank: " << my_rank << " increase tasks" << std::endl;

		//std::cout << "current_G edges: " << current_G->GetNumEdges() << " vertices "<< current_G->GetNumVertices() << std::endl;
		//std::cout << "u: " << u << " v: " << v << std::endl;

		std::vector<Branch> new_branches(2);

		//std::cout << " ub_og " << std::endl;
		int lb_og = _clique_strat.FindClique(*current_G);
		unsigned short ub_og;
		std::cout << "lb_og " << lb_og << " ub_og " << ub_og << std::endl;
		_color_strat.Color(*current_G, ub_og);
		//std::cout << "lb_og (again)" << lb_og << " ub_og (again)" << ub_og << std::endl;
		Log_par("[Branch 0] Original "
			"lb = " + std::to_string(lb_og) +
			", ub = " + std::to_string(ub_og),
			depth);
		
		// MergeVertices
		auto G1 = current_G->Clone();
		//std::cout << "Rank: " << my_rank << " copy graph" << std::endl;
		G1->MergeVertices(u, v);
		//std::cout << "Rank: " << my_rank << " merge vertices" << std::endl;
		int lb1 = _clique_strat.FindClique(*G1);
		unsigned short ub1;
		_color_strat.Color(*G1, ub1);
		Log_par("[Branch 1] (Merge u, v) "
			"lb = " + std::to_string(lb1) +
			", ub = " + std::to_string(ub1),
			depth);

		if (lb1 < best_ub.load()) {  // check to avoid cuncurrent access
			new_branches[0] = Branch(std::move(G1), lb1, ub1, 1);
		}
		// AddEdge
		auto G2 = current_G->Clone();
		G2->AddEdge(u, v);
		int lb2 = _clique_strat.FindClique(*G2);
		unsigned short ub2;
		_color_strat.Color(*G2, ub2);
		Log_par("[Branch 2] (Add edge u-v) "
			"lb = " + std::to_string(lb2) +
			", ub = " + std::to_string(ub2),
			depth);

		if ((lb2 < best_ub.load()) &&
			(lb2 < new_branches[0].ub)) {  // check to avoid cuncurrent access
			new_branches[1] = Branch(std::move(G2), lb2, ub2, 1);
		}

		active_tasks.fetch_sub(1);
		std::cout << "Rank: " << my_rank << " decrease tasks" << std::endl;
	
		{
			std::lock_guard<std::mutex> lock(queue_mutex);
			for (auto& branch : new_branches) {
				if (branch.g) queue.push(std::move(branch));
			}
		}

		// Update local sbest_ub
		unsigned short previous_best_ub = best_ub.load();
		best_ub.store(std::min({previous_best_ub, new_branches[0].ub, new_branches[1].ub}));
		Log_par("[UPDATE] Updated best_ub: " + std::to_string(best_ub), depth);

	active_tasks.fetch_add(1);
	*/
}

/*CHANGE:   now ub and best_ub are atomic variables(avoids concurrent access)
			added mutex to avoid concurrent access to the queue
			implemented system to send work to workers(line 421)
*/
// TODO: Add a flag to terminate the loop and avoid infinite loop
int BranchNBoundPar::Solve(Graph& g, int timeout_seconds, int iteration_threshold) {
	// Intitialize MPI
	//MPI_Init(NULL, NULL);

	int my_rank;
	int p;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Start the timeout timer
	auto global_start_time = MPI_Wtime();
	// Initialize big enough best_ub for all processes.
	std::atomic<unsigned short> best_ub = USHRT_MAX;
	std::mutex queue_mutex;	 // avoid concurrent access to the queue
	BranchQueue queue;

	// TODO: possibility to reduce initial time by distributing
	// work(generate 2 branches, keep one send the other, i think this one
	// is more scalable) to workers or using openMP
	if (my_rank == 0) {
		// Initialize bounds
		int lb = _clique_strat.FindClique(g);
		unsigned short ub;
		_color_strat.Color(g, ub);
		best_ub = ub;

		// Log initial bounds
		Log("[INIT] Initial bounds: lb = " + std::to_string(lb) +
		    ", ub = " + std::to_string(ub));

		queue.push(Branch(g.Clone(), lb, ub, 1));	// Initial branch with depth 1

		// Start computation tree
		Log("========== START ==========", 0);
		Log("Root: (lb = " + std::to_string(lb) +
			", ub = " + std::to_string(ub) + ")",
		    0);

		// Get enough branches to distribute amongst workers
		while (queue.size() < p - 1) {
			// Dequeue the next branch
			Branch current = std::move(const_cast<Branch&>(queue.top()));
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
			auto [u, v] = _branching_strat.ChooseVertices(*current_G);
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
			best_ub = std::min({best_ub.load(), ub1, ub2});
			Log("[UPDATE] Updated best_ub: " +
				std::to_string(best_ub),
			    current.depth);
		}
		// assume we know the number of workers (p) and we have p-1
		// branches, we send to each worker a branch

		// std::cout << "Queue size" << queue.size() << std::endl;
		for (int i = 1; i <= p - 1; i++) {
			Branch branch_to_send = std::move(const_cast<Branch&>(queue.top()));
			queue.pop();
			// send branch to worker i
			sendBranch(branch_to_send, i, TAG_WORK, MPI_COMM_WORLD);
		}
		// std::cout << "Queue size" << queue.size() << std::endl;

		Log("[PARALLELISATION START]");
	} else {
		MPI_Status status_recv;
		Branch branch_recv;

		// recv work
		branch_recv = recvBranch(0, TAG_WORK, MPI_COMM_WORLD);

		/*
		Dimacs dimacs;
		std::string file_name = "10_vertices_graph.col";

		if (!dimacs.load(file_name.c_str())) {
			std::cout << dimacs.getError() << std::endl;
			return 1;
		}

		CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);

		if(graph->isEqual(*branch_recv.g)){
			std::cout << "Graphs are equal" << std::endl;
		}
		else{
			std::cout << "Graphs are not equal" << std::endl;
		}

		*/

		std::atomic<unsigned short> best_ub = branch_recv.ub;

		//std::cout << "ub_recv: " << _clique_strat.FindClique(*branch_recv.g) << std::endl;

		//std::cout << "Worker: " << my_rank << std::endl;

		queue.push(std::move(branch_recv));

		//std::cout << "Pushed to queue" << std::endl;
	}

	// OpenMP Parallel Region
	/*
	idea is assign specific threads to specific tasks, in particular the
	first three threads are used for gathering best_ub, checking for
	termination and listening for work requests, then the remaining threads are used
	for -one to keep popping from the queue, work_stealing, generate
	omp_tasks(method create_task) and add new branches to the queue -others
	to compute the omp_tasks
	*/
	omp_set_num_threads(5);
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();

		// Both master's and worker's thread 0 goes in here.
		if (tid == 0) {
			//it must be thread 0 otherwise it will not work
			// Checks if solution has been found or timeout. 
			thread_0_terminator(my_rank, p, global_start_time, timeout_seconds);
			std::cout << "Rank: " << my_rank << " Exited thread_0_terminator func." << std::endl;
		}
		// Both master's and worker's thread 1 goes in here.
		if (tid == 1) {
			// Updates (gathers) best_ub from time to time.
			thread_1_solution_gatherer(p, best_ub, my_rank);
			std::cout << "Rank: " << my_rank << " Exited thread_1_solution_gatherer func." << std::endl;
		}
		// Only worker processes go in here.
		if (my_rank != 0) {
			// worker thread used for listening for requests about
			// work stealing
			if (tid == 2) {
				thread_2_listen_for_requests(queue_mutex, queue, my_rank);
				std::cout << "Rank: " << my_rank << " Exited thread_2_listen_for_requests func." << std::endl;
			}

			std::unique_ptr<Graph> current_G;
			int current_lb;
			unsigned short current_ub;
			unsigned int u, v;

				//std::cout << "Rank: " << my_rank << " Starting work" << std::endl;
				Branch current;

				while (!terminate_flag.load(std::memory_order_relaxed)) {
					bool has_work = false;

					{
						//std::cout << "Rank: " << my_rank << " checking queue" << std::endl;
						std::lock_guard<std::mutex> lock(queue_mutex);
						if (!queue.empty()) { 
							current = std::move(const_cast<Branch&>(queue.top()));
							queue.pop();
							has_work = true;
						}
					}

					//std::cout << "Rank: " << my_rank << " queue checked" << std::endl;

					if (!has_work) {
						//std::cout << "Rank: " << my_rank << " no work in queue" << std::endl;
						if(!work_stealing(my_rank, p, queue, queue_mutex, current)){
							continue;
						}
					}

					current_G = std::move(current.g);
					current_lb = current.lb;
					current_ub = current.ub;

					Log_par("Processing node: lb = " + std::to_string(current_lb) +
						    ", ub = " + std::to_string(current_ub), current.depth);

					// find solution
					if (current_lb == current_ub) {
						MPI_Send(&current_ub, 1, MPI_INT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD);  // check if it is correct
						Log_par(
						    "[FOUND] Chromatic number "
						    "found: " + std::to_string(current_lb),current.depth);
						Log_par("========== END ==========", 0);
						continue;
					}

					// Prune
					if (current_lb >= best_ub.load()) {
						Log_par(
						    "[PRUNE] Branch pruned at "
						    "depth " + std::to_string(current.depth) +
							": lb = " + std::to_string(current_lb) +
							" >= best_ub = " + std::to_string(best_ub.load()),
						    current.depth);
						continue;
					}

					//std::cout << "Rank: " << my_rank << " starting branching" << std::endl;

					auto [u, v] = _branching_strat.ChooseVertices(*current_G);
					Log_par("Branching on vertices: u = " + std::to_string(u) +
							", v = " + std::to_string(v),
					    	current.depth, true);

					//std::cout << "Rank: " << my_rank << " choose " << u << " " << v << std::endl;

					if (u == -1 || v == -1) {
						unsigned short chromatic_number = std::min<unsigned short>(current_G->GetNumVertices(), best_ub.load()); //check this
						MPI_Send(&chromatic_number, 1, MPI_UNSIGNED_SHORT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD);  // check if it is correct
						Log_par("Graph is complete. "
						    	"Chromatic number = " + std::to_string(chromatic_number),
						    	current.depth);
						Log_par("========== END ==========", 0);
						continue;
					}

					// TODO: check for other approaches to
					// avoid generating too much tasks(this
					// will wait untill all child tasks are
					// completed) avoid generating too much
					// tasks

					//std::cout << "Rank: " << my_rank << " checking # tasks" << std::endl;

					//std::cout << "Rank: " << my_rank << " generating tasks" << std::endl;
					// generate tasks and update queue
					int current_depth = current.depth;

					//int lb_og_pre = _clique_strat.FindClique(*current_G);
					//unsigned short ub_og_pre;
					//_color_strat.Color(*current_G, ub_og_pre);
					//std::cout << "lb_og_pre " << lb_og_pre << " ub_og_pre " << ub_og_pre << std::endl;

					//auto local_g = current_G->Clone();

					//std::cout << "Rank: " << my_rank << " start task" << std::endl;
					//active_tasks.fetch_add(1);
					//std::cout << "Rank: " << my_rank << " increase tasks" << std::endl;

					//std::cout << "current_G edges: " << current_G->GetNumEdges() << " vertices "<< current_G->GetNumVertices() << std::endl;
					//std::cout << "u: " << u << " v: " << v << std::endl;

					//std::vector<Branch> new_branches(2);

					//std::cout << " ub_og " << std::endl;
					/*
					int lb_og = _clique_strat.FindClique(*current_G);
					unsigned short ub_og;
					//std::cout << "lb_og " << lb_og << " ub_og " << ub_og << std::endl;
					_color_strat.Color(*current_G, ub_og);
					//std::cout << "lb_og (again)" << lb_og << " ub_og (again)" << ub_og << std::endl;
					Log_par("[Branch 0] Original "
						"lb = " + std::to_string(lb_og) +
						", ub = " + std::to_string(ub_og),
						current_depth);
					*/
					
					// MergeVertices
					auto G1 = current_G->Clone();
					//std::cout << "Rank: " << my_rank << " copy graph" << std::endl;
					G1->MergeVertices(u, v);
					//std::cout << "Rank: " << my_rank << " merge vertices" << std::endl;
					int lb1 = _clique_strat.FindClique(*G1);
					unsigned short ub1;
					_color_strat.Color(*G1, ub1);
					Log_par("[Branch 1] (Merge u, v) "
						"lb = " + std::to_string(lb1) +
						", ub = " + std::to_string(ub1),
						current_depth);

					if (lb1 < best_ub.load()) { // check to avoid cuncurrent access
						std::lock_guard<std::mutex> lock(queue_mutex);
						queue.push(Branch(std::move(G1), lb1, ub1, current_depth + 1));
					}
					// AddEdge
					auto G2 = current_G->Clone();
					G2->AddEdge(u, v);
					int lb2 = _clique_strat.FindClique(*G2);
					unsigned short ub2;
					_color_strat.Color(*G2, ub2);
					Log_par("[Branch 2] (Add edge u-v) "
						"lb = " + std::to_string(lb2) +
						", ub = " + std::to_string(ub2),
						current_depth);

					if (lb2 < best_ub.load() && (lb2 < ub1)) { 
						std::lock_guard<std::mutex> lock(queue_mutex);
						queue.push(Branch(std::move(G2), lb2, ub2, current_depth + 1));
					}

					//active_tasks.fetch_sub(1);
					//std::cout << "Rank: " << my_rank << " decrease tasks" << std::endl;

					// Update local sbest_ub
					unsigned short previous_best_ub = best_ub.load();
					best_ub.store(std::min({previous_best_ub, ub1, ub2}));
					Log_par("[UPDATE] Updated best_ub: " + std::to_string(best_ub.load()), current_depth);

					
				//active_tasks.fetch_add(1);
				}
			}
		//# pragma omp barrier
	}
	std::cout << "Rank: " << my_rank << " Finalizing" << std::endl;
	// End execution
	return best_ub;
}

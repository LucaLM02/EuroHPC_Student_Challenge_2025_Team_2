#include "branch_n_bound_par.hpp"

#define ALLGATHER_WAIT_TIME 5	   // Sleep time for MPI_Allgather
#define TIMEOUT_CHECK_WAIT_TIME 5  // Sleep time for timeout checker

// tags for MPI communication
#define TAG_WORK_REQUEST 1
#define TAG_WORK_RESPONSE 2
#define TAG_WORK 3
#define TAG_SOLUTION_FOUND 4
#define TAG_IDLE 5

using BranchQueue = std::priority_queue<Branch, std::vector<Branch>>;
std::atomic<bool> terminate_flag = false;
std::mutex queue_mutex;	 // avoid concurrent access to the queue
std::mutex log_mutex;

std::mutex branching_mutex;
std::mutex task_mutex;

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


void BranchNBoundPar::Log_par(const std::string& message, int depth = 0, bool is_branching = false) {
	std::lock_guard<std::mutex> lock(log_mutex);
	if (_log_file_branches.is_open()) {
		int rank = 0, thread_id = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		thread_id = omp_get_thread_num();

		// Indentation based on depth
		std::string indentation(depth * 2, ' ');
		// TODO: Logging everything with this critical section is not efficient
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
	std::vector<char> buffer = b.serialize();
	int size = buffer.size();

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
	int size = 0;

	MPI_Recv(&size, 1, MPI_INT, source, tag, comm, &status);

	std::vector<char> buffer(size, 0);
	MPI_Recv(buffer.data(), size, MPI_BYTE, source, tag, comm, &status);
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
void thread_0_terminator(int my_rank, int p, int global_start_time, int timeout_seconds, double &optimum_time) {
	int solution_found = 0;
	int timeout_signal = 0;

	std::vector<int> idle_status(p, 0); // Array to keep track of idle status of workers
	while (true) {
		if (my_rank == 0) {
			// Master listens for solution found (Non-blocking)
			MPI_Status status;
			int flag = 0;
			MPI_Iprobe(MPI_ANY_SOURCE, TAG_SOLUTION_FOUND, MPI_COMM_WORLD, &flag, &status);
			// Check if a solution is being communicated
			if (flag) {
				unsigned short solution = 0;
				MPI_Recv(&solution, 1, MPI_UNSIGNED_SHORT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				solution_found = 1;
			}
			// Check if timeout is reached, broadcast timeout signal
			if (MPI_Wtime() - global_start_time >= timeout_seconds) {
				timeout_signal = 1;
			}

			// Listen for idle status updates from workers
            MPI_Iprobe(MPI_ANY_SOURCE, TAG_IDLE, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int worker_idle_status = 0;
                MPI_Recv(&worker_idle_status, 1, MPI_INT, status.MPI_SOURCE, TAG_WORK_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                idle_status[status.MPI_SOURCE] = worker_idle_status;
            }

            // Check if all workers are idle
            if (std::all_of(idle_status.begin(), idle_status.end(), [](int status) { return status == 1; })) {
                solution_found = 1;
				optimum_time = MPI_Wtime() - global_start_time;
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
void thread_1_solution_gatherer(int p, std::atomic<unsigned short>& best_ub, int& my_rank) { 
    std::vector<unsigned short> all_best_ub(p);
    auto last_gather_time = MPI_Wtime();
    MPI_Request request;
	
    while (!terminate_flag.load(std::memory_order_relaxed)) {
        auto current_time = MPI_Wtime();
        auto elapsed_time = current_time - last_gather_time;

        if (elapsed_time >= ALLGATHER_WAIT_TIME) {
            unsigned short local_best_ub = best_ub.load(); // safe read

            // Start non-blocking allgather
            MPI_Iallgather(&local_best_ub, 1, MPI_UNSIGNED_SHORT, all_best_ub.data(), 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD, &request);

            // Wait for completion with timeout handling (or simply test it periodically)
            MPI_Status status;
            while (true) {
                int flag = 0;
                MPI_Test(&request, &flag, &status);
                if (flag) break;  // The operation is completed

                // If termination flag is set, cancel the request to avoid deadlock
                if (terminate_flag.load(std::memory_order_relaxed)) {
                    return;
                }

                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }

            // Update the best upper bound for other threads in shared memory
            best_ub.store(*std::min_element(all_best_ub.begin(), all_best_ub.end()));  // safe write

            // Reset the timer
            last_gather_time = current_time;
        }
        
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}


/**
 * thread_2_employer - Listens for work requests from other worker processes and
 * distributes available work from the local queue. This function is called by
 * the third thread of the worker processes.
 *
 * This function is responsible for handling work requests from other worker
 * processes. It listens for incoming work requests using non-blocking MPI
 * communication. If a work request is received and there is work available in
 * the local queue, it sends a branch to the requesting process. If no work is
 * available, it sends a response indicating that no work is available.
 *
 * Parameters:
 *   queue_mutex (std::mutex&) : Mutex to protect concurrent access to the work
 * queue.
 *   queue (BranchQueue&)      : The local work queue containing branches to be
 * processed.
 *   my_rank (int&)            : The rank of the current process.
 */
void thread_2_employer(std::mutex& queue_mutex, BranchQueue& queue, int& my_rank) {
    MPI_Status status;
    int request_signal = 0;
    while (!terminate_flag.load(std::memory_order_relaxed)) {
        // Listen for a request for work from other workers.
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_WORK_REQUEST, MPI_COMM_WORLD, &request_signal, &status);
        int destination_rank = status.MPI_SOURCE;
        int response = 0;
        // Check if a solution is being communicated
        if (request_signal) {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (!queue.empty()) {
                response = 1;
                MPI_Send(&response, 1, MPI_INT, destination_rank, TAG_WORK_RESPONSE, MPI_COMM_WORLD);
                Branch branch = std::move(const_cast<Branch&>(queue.top()));
                queue.pop();
                lock.unlock();
                sendBranch(branch, destination_rank, TAG_WORK, MPI_COMM_WORLD);
            } else {
                lock.unlock();
                MPI_Send(&response, 1, MPI_INT, destination_rank, TAG_WORK_RESPONSE, MPI_COMM_WORLD);
            }
        }
    }
}

/**
 * request_work - Requests work from other worker processes when the local queue is empty.
 *
 * This function is responsible for requesting work from other worker processes
 * when the local queue is empty. It uses non-blocking MPI communication to send
 * a work request to a randomly selected worker and waits for a response. If work
 * is available, it receives the branch and adds it to the local queue.
 *
 * Parameters:
 *   my_rank (int)            : The rank of the current process.
 *   p (int)                  : The total number of processes in the MPI communicator.
 *   queue (BranchQueue&)     : The local work queue containing branches to be processed.
 *   queue_mutex (std::mutex&): Mutex to protect concurrent access to the work queue.
 *   current (Branch&)        : The branch object to store the received work.
 *
 * Returns:
 *   bool : True if work was successfully received and added to the queue, false otherwise.
 */
bool request_work(int my_rank, int p, BranchQueue& queue, std::mutex& queue_mutex, Branch& current) {
    int target_worker = (rand() % (p - 1)) + 1; // Randomly select a worker to request work from
    MPI_Status status;
    int response = 0;
    MPI_Request request;

    MPI_Isend(nullptr, 0, MPI_INT, target_worker, TAG_WORK_REQUEST, MPI_COMM_WORLD, &request);
    MPI_Irecv(&response, 1, MPI_INT, target_worker, TAG_WORK_RESPONSE, MPI_COMM_WORLD, &request);

    while (true) {
        int flag = 0;
        MPI_Test(&request, &flag, &status);
        if (flag) break;  // The operation is completed

        // If termination flag is set, cancel the request to avoid deadlock
        if (terminate_flag.load(std::memory_order_relaxed)) {
            MPI_Cancel(&request);
            return false;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    if (response == 1) { // Work is available
        current = recvBranch(target_worker, TAG_WORK, MPI_COMM_WORLD);
        std::lock_guard<std::mutex> lock(queue_mutex);
        queue.push(std::move(current));
        return true;
    }
    return false;
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
 /*
void BranchNBoundPar::create_task(
    std::unique_ptr<Graph> current_G, int u, int v,
    CliqueStrategy& _clique_strat, ColorStrategy& _color_strat,
	std::atomic<unsigned short> & best_ub, int const & depth, int my_rank) {

		//active_tasks.fetch_add(1);
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

				if (lb1 < best_ub.load()) { // check to avoid cuncurrent access
					std::lock_guard<std::mutex> lock(queue_mutex);
					queue.push(Branch(std::move(G1), lb1, ub1, depth + 1));
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

				if (lb2 < best_ub.load() && (lb2 < ub1)) { 
					std::lock_guard<std::mutex> lock(queue_mutex);
					queue.push(Branch(std::move(G2), lb2, ub2, depth + 1));
				}

				//active_tasks.fetch_sub(1);
				//std::cout << "Rank: " << my_rank << " decrease tasks" << std::endl;

				// Update local sbest_ub
				unsigned short previous_best_ub = best_ub.load();
				best_ub.store(std::min({previous_best_ub, ub1, ub2}));
				Log_par("[UPDATE] Updated best_ub: " + std::to_string(best_ub.load()), depth);

		//active_tasks.fetch_sub(1);
}
*/

int BranchNBoundPar::Solve(Graph& g, double &optimum_time, int timeout_seconds) {
	optimum_time  = -1.0;

	BranchQueue queue;
	int my_rank;
	int p;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Start the timeout timer.
	auto global_start_time = MPI_Wtime();
	// Initialize big enough best_ub for all processes.
	std::atomic<unsigned short> best_ub = USHRT_MAX;

	MPI_Status status_recv;
	Branch branch_recv;

	// Receive work from previous worker if not rank 0 (starting rank), to start working.
	if (my_rank>0) {
		branch_recv = recvBranch(my_rank-1, TAG_WORK, MPI_COMM_WORLD);
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

	queue.push(std::move(branch_recv));

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
	omp_set_num_threads(4);
	#pragma omp parallel default(shared)
	{
		int tid = omp_get_thread_num();

		if (tid == 0) { // Checks if solution has been found or timeout. 
			thread_0_terminator(my_rank, p, global_start_time, timeout_seconds, optimum_time);
			std::cout << "Rank: " << my_rank << " Exited thread_0_terminator func." << std::endl;
		}else if (tid == 1) { // Updates (gathers) best_ub from time to time.
			thread_1_solution_gatherer(p, best_ub, my_rank);
			std::cout << "Rank: " << my_rank << " Exited thread_1_solution_gatherer func." << std::endl;
		}else if (tid == 2) { // Employer thread employs workers by answering their work requests
			thread_2_employer(queue_mutex, queue, my_rank);
		}else if (tid == 3) { // TODO: Let more threads do these computations in parallel
			Branch current;
			NeighboursBranchingStrategy* branching_strategy_local = dynamic_cast<NeighboursBranchingStrategy*>(&_branching_strat);
			FastCliqueStrategy* clique_strategy_local = dynamic_cast<FastCliqueStrategy*>(&_clique_strat);
			GreedyColorStrategy* color_strategy_local = dynamic_cast<GreedyColorStrategy*>(&_color_strat);

				std::cout << "Rank: " << my_rank << " Starting... " << std::endl;

				// If rank 0, initialize first branch.
				if (my_rank == 0) {
					// Initialize bounds
					int lb = _clique_strat.FindClique(g);
					unsigned short ub;
					_color_strat.Color(g, ub);
					best_ub = ub;
			
					// Log initial bounds
					Log_par("[INIT] Initial bounds: lb = " + std::to_string(lb) +
						", ub = " + std::to_string(ub));

					std::lock_guard<std::mutex> lock(queue_mutex);
					queue.push(Branch(g.Clone(), lb, ub, 1));	// Initial branch with depth 1
				}

				bool distributed_work = false; // Signals when work distribution phase ends.

				while (!terminate_flag.load()) {
					bool has_work = false;
					{
						std::lock_guard<std::mutex> lock(queue_mutex);
						if (!queue.empty()) { 
							current = std::move(const_cast<Branch&>(queue.top()));
							queue.pop();
							has_work = true;
						}
					}

					// If no work and already passed the initial distributing phase, request work.
					if (!has_work && distributed_work) {
						#pragma omp single // Only a single thread asks for work.
						{
						// Notify the root process that this worker is idle
						int idle_status = 1;
						MPI_Send(&idle_status, 1, MPI_INT, 0, TAG_IDLE, MPI_COMM_WORLD);
						// Start requesting work.
						std::cout << "Rank: " << my_rank << " requesting work..." << std::endl;
						while (!request_work(my_rank, p, queue, queue_mutex, current)) {
							if (terminate_flag.load()) break;
							std::this_thread::sleep_for(std::chrono::milliseconds(10));
						}
						// Work received. Notify the root process that this worker is not idle anymore.
						idle_status = 0;
						MPI_Send(&idle_status, 1, MPI_INT, 0, TAG_IDLE, MPI_COMM_WORLD);				
						}
						continue;
					}

					auto current_G = std::move(current.g);
					int current_lb = current.lb;
					unsigned short current_ub = current.ub;

					Log_par("Processing node: lb = " + std::to_string(current_lb) +
							", ub = " + std::to_string(current_ub), current.depth);

					// found solution
					if (current_lb == current_ub) {
						/*
						MPI_Send(&current_ub, 1, MPI_UNSIGNED_SHORT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD);  // check if it is correct
						Log_par(
							"[FOUND] Chromatic number "
							"found: " + std::to_string(current_lb),current.depth);
						Log_par("========== END ==========", 0);
						*/
						best_ub.store(current_ub);
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

					// Start branching 
					std::unique_lock<std::mutex> lock_branching(branching_mutex);
					auto [u, v] = _branching_strat.ChooseVertices(*current_G);
					lock_branching.unlock();
					Log_par("Branching on vertices: u = " + std::to_string(u) +
							", v = " + std::to_string(v),
							current.depth, true);

					if (u == -1 || v == -1) {
						best_ub.store(std::min<unsigned short>(current_G->GetNumVertices(), best_ub.load()));
						/*
						MPI_Send(&chromatic_number, 1, MPI_UNSIGNED_SHORT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD);  // check if it is correct
						Log_par("Graph is complete. "
								"Chromatic number = " + std::to_string(chromatic_number),
								current.depth);
						Log_par("========== END ==========", 0);
						*/
						continue;
					}
					// generate tasks and update queue
					std::unique_lock<std::mutex> lock_task(task_mutex);
					auto G1 = current_G->Clone();
					G1->MergeVertices(u, v);
					int lb1 = _clique_strat.FindClique(*G1);
					unsigned short ub1;
					_color_strat.Color(*G1, ub1);
					Log_par("[Branch 1] (Merge u, v) "
							"lb = " + std::to_string(lb1) +
							", ub = " + std::to_string(ub1),
							current.depth);

					if (lb1 < best_ub.load()) { // check to avoid cuncurrent access
						std::lock_guard<std::mutex> lock(queue_mutex);
						queue.push(Branch(std::move(G1), lb1, ub1, current.depth + 1));
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
							current.depth);

					if (lb2 < best_ub.load() && (lb2 < ub1)) { 
						std::lock_guard<std::mutex> lock(queue_mutex);
						queue.push(Branch(std::move(G2), lb2, ub2, current.depth + 1));
					}
					if (!distributed_work && my_rank<p-1)
					{
						std::cout << "Rank: " << my_rank << " Distributing work" << std::endl;
						std::lock_guard<std::mutex> lock(queue_mutex);{
						current = std::move(const_cast<Branch&>(queue.top()));
						queue.pop();}
						sendBranch(current, my_rank+1, TAG_WORK, MPI_COMM_WORLD);
						distributed_work = true;
					}
					lock_task.unlock();
					// Update local sbest_ub
					unsigned short previous_best_ub = best_ub.load();
					best_ub.store(std::min({previous_best_ub, ub1, ub2}));
					Log_par("[UPDATE] Updated best_ub: " + std::to_string(best_ub.load()), current.depth);
			}
		}
	}
		std::cout << "Rank: " << my_rank << " Finalizing. " << std::endl;
		// End execution
		return best_ub;
	}
		














				/*
			#pragma omp single nowait
			{

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
						/*
						MPI_Send(&current_ub, 1, MPI_UNSIGNED_SHORT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD);  // check if it is correct
						Log_par(
						    "[FOUND] Chromatic number "
						    "found: " + std::to_string(current_lb),current.depth);
						Log_par("========== END ==========", 0);
						*/
						/*
						best_ub.store(current_ub);
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
						best_ub.store(std::min<unsigned short>(current_G->GetNumVertices(), best_ub.load()));
						/*
						MPI_Send(&chromatic_number, 1, MPI_UNSIGNED_SHORT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD);  // check if it is correct
						Log_par("Graph is complete. "
						    	"Chromatic number = " + std::to_string(chromatic_number),
						    	current.depth);
						Log_par("========== END ==========", 0);
						*/
						//continue;
					//}

					// TODO: check for other approaches to
					// avoid generating too much tasks(this
					// will wait untill all child tasks are
					// completed) avoid generating too much
					// tasks

					//std::cout << "Rank: " << my_rank << " checking # tasks" << std::endl;
					/*
					if (active_tasks.load() > max_tasks) {
						#pragma omp taskwait
					}

					//std::cout << "Rank: " << my_rank << " generating tasks" << std::endl;
					// generate tasks and update queue
					int current_depth = current.depth;
					*/
					/*
					int lb_og_pre = _clique_strat.FindClique(*current_G);
					unsigned short ub_og_pre;
					_color_strat.Color(*current_G, ub_og_pre);
					//std::cout << "lb_og_pre " << lb_og_pre << " ub_og_pre " << ub_og_pre << std::endl;

					auto local_g = current_G->Clone();
					*/

					//Graph* local_g = current_G->Clone().release();
					/*
					std::unique_lock<std::mutex> lock_current_G(current_G_mutex);
					std::cout << "num vertices " << current_G->GetNumVertices() << " edges " <<  current_G->GetNumEdges() << " u " << u << " v " << v << std::endl;
					graph_queue.emplace(std::move(current_G), u, v);
					lock_current_G.unlock();

					#pragma omp task default(shared)	
					{
						std::unique_lock<std::mutex> lock_current_G(current_G_mutex);
						GraphEntry local = std::move(graph_queue.front());
						std::unique_ptr<Graph> local_g = std::move(local.graph);
						int u_local = local.u;
						int v_local = local.v;
						graph_queue.pop();
						std::cout << "num vertices " << local_g->GetNumVertices() << " edges " <<  local_g->GetNumEdges() << " u " << u_local << " v " << v_local << std::endl;
						lock_current_G.unlock();
						//ptr.reset(local_g);
						create_task(std::move(local_g), u_local, v_local, _clique_strat, _color_strat, best_ub, current_depth, my_rank);
					}

				}
			}
				
		}

		//# pragma omp barrier
	}
	std::cout << "Rank: " << my_rank << " Finalizing" << std::endl;
	// End execution
	return best_ub;
}
			*/
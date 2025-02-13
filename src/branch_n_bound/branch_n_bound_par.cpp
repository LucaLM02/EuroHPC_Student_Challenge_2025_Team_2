#include "branch_n_bound_par.hpp"

#define ALLGATHER_WAIT_TIME 10	   // Sleep time for MPI_Allgather
#define TIMEOUT_CHECK_WAIT_TIME 5  // Sleep time for timeout checker

//tags for MPI communication
#define TAG_WORK_REQUEST 1
#define TAG_WORK_RESPONSE 2
#define TAG_WORK 3
#define TAG_SOLUTION_FOUND 4

using BranchQueue = std::priority_queue<Branch, std::vector<Branch>>;

/*
* FOR JAN
* TODO == suggested changes/improvements
* CHANGE == changes made
*/

//CHANGE: fixed with casting start time as double
//TODO: check if its better using only mpi_wtime instead of chrono
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

//CHANGE: added thread_id to log and make it thread safe(by using critical section)
void BranchNBoundPar::Log_par(const std::string& message, int depth = 0,
	bool is_branching = false) {
	if (_log_file.is_open()) {
		int rank=0, thread_id=0;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		thread_id = omp_get_thread_num();

		// Indentation based on depth
		std::string indentation(depth * 2, ' ');

		#pragma omp critical
		{
			if (is_branching) {
				_log_file << indentation << "[Rank " << rank 
							<< " | Thread " << thread_id << "] "
							<< "[Depth " << depth << "] Branching on u = " 
							<< message << std::endl;
			} else {
				_log_file << indentation << "[Rank " << rank 
							<< " | Thread " << thread_id << "] "
							<< message << std::endl;
			}
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

 //CHANGE: all_best_ub is now a vector of unsigned short, best_ub is now an atomic variable(helps to avoid cuncurrent access)
void thread_0_solution_gatherer(int p, std::atomic<unsigned short>& best_ub) {
	std::vector<unsigned short> all_best_ub(p);
	// Start timer for best_ub exchange interval.
	time_t start_time = time(NULL);

	while (1) { //TODO: use flag to avoid infinite loop

		// Wait for the time threshold. Allgather in done every
		// ALLGATHER_WAIT_TIME seconds.
		sleep(ALLGATHER_WAIT_TIME);
		unsigned short local_best_ub = best_ub.load(); // safe read

		// Gather best_ub from other workers
		MPI_Allgather(&local_best_ub, 1, MPI_UNSIGNED_SHORT, all_best_ub.data(), 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);

		// Update the best upper bound for other threads in
		// shared memory
		best_ub.store(*std::min_element(all_best_ub.begin(), all_best_ub.end())); //safe write

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
//TODO: evaluate using iprobe instead of recv to reduce comunication overhead
//CHANGE: added mutex to avoid concurrent access to the queue, added response to avoid sending work to workers if there is no work available(avoid also deadlock in workers)
void thread_2_listen_for_requests(std::mutex &queue_mutex, BranchQueue &queue) {
	MPI_Status status;
	int request_signal;
	while (true) {
        // Listen for a request for work from other workers.
        MPI_Recv(&request_signal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_WORK_REQUEST, MPI_COMM_WORLD, &status);

		int destination_rank = status.MPI_SOURCE;
		int response = 0; 

		std::lock_guard<std::mutex> lock(queue_mutex);
        if (!queue.empty()) {
			response = 1;
			MPI_Send(&response, 1, MPI_INT, destination_rank, TAG_WORK_REQUEST, MPI_COMM_WORLD);

            Branch branch = std::move(const_cast<Branch&>(queue.top()));
            queue.pop();

            MPI_Send(&branch, sizeof(Branch), MPI_BYTE, destination_rank, TAG_WORK, MPI_COMM_WORLD);
        }else{
			MPI_Send(&response, 1, MPI_INT, destination_rank, TAG_WORK_REQUEST, MPI_COMM_WORLD);
		}
	}
}

/**
 * Creates a task for the OpenMP parallel region to execute.
 *  
 * 
 * @param active_tasks The number of active tasks in the parallel region.
 * @param current_G The current graph.
 * @param u The first vertex to merge or add an edge.
 * @param v The second vertex to merge or add an edge.
 * @param _clique_strat The clique strategy.
 * @param _color_strat The color strategy.
 * @param new_branches The new branches to add to the queue.
 * @param task_type The type of task to create (1 = MergeVertices, 2 = AddEdge).
 * @param best_ub The best upper bound.
 */
void BranchNBoundPar::create_task(std::atomic<int>& active_tasks, Graph* current_G, int u, int v,
CliqueStrategy& _clique_strat, ColorStrategy& _color_strat, 
std::vector<Branch>& new_branches, int task_type, std::atomic<unsigned short>const &best_ub, int const &depth) {

	#pragma omp atomic
    active_tasks++;

	if (task_type == 1) {
        // MergeVertices
        auto G1 = current_G->Clone();
        G1->MergeVertices(u, v);
        int lb1 = _clique_strat.FindClique(*G1);
        unsigned short ub1;
        _color_strat.Color(*G1, ub1);
		Log_par("[Branch 1] (Merge u, v) lb = " + std::to_string(lb1) +
			", ub = " + std::to_string(ub1),
		    depth);

        if (lb1 < best_ub.load()) { //check to avoid cuncurrent access
            new_branches[0] = Branch(std::move(G1), lb1, ub1, 1); 
        }
    } else if (task_type == 2) {
        // AddEdge
        auto G2 = current_G->Clone();
        G2->AddEdge(u, v);
        int lb2 = _clique_strat.FindClique(*G2);
        unsigned short ub2;
        _color_strat.Color(*G2, ub2);
		Log_par("[Branch 2] (Add edge u-v) lb = " + std::to_string(lb2) +
			", ub = " + std::to_string(ub2),
		    depth);

        if ((lb2 < best_ub.load()) && (lb2 < new_branches[0].ub)) { //check to avoid cuncurrent access
            new_branches[1] = Branch(std::move(G2), lb2, ub2, 1);
        }
	}

    #pragma omp atomic
    active_tasks--;
}

/*CHANGE:   now ub and best_ub are atomic variables(avoids concurrent access)
			added mutex to avoid concurrent access to the queue
			implemented system to send work to workers(line 421)
*/
// TODO: Add a flag to terminate the loop and avoid infinite loop
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
	std::atomic<unsigned short> best_ub = USHRT_MAX;

	std::mutex queue_mutex; //avoid concurrent access to the queue
	BranchQueue queue;

	std::atomic<int> active_tasks = 0;
	int max_tasks = omp_get_max_threads()-3; //limit number of tasks (maybe we can increase this number)

	//TODO: possibility to reduce initial time by distributing work(generate 2 branches, keep one send the other, i think this one is more scalable) to workers 
	//or using openMP
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
			best_ub = std::min({best_ub.load(), ub1, ub2});
			Log("[UPDATE] Updated best_ub: " +
				std::to_string(best_ub),
			    current.depth);
		}

		//assume we know the number of workers (p) and we have p-1 branches, we send to each worker a branch
		for (int i = 1; i <= p; i++) {
			Branch branch_to_send = std::move(const_cast<Branch&>(queue.top()));
			queue.pop();
	
			// send branch to worker i
			MPI_Send(&branch_to_send, sizeof(Branch), MPI_BYTE, i, TAG_WORK, MPI_COMM_WORLD);
		}
		
		Log("[PARALLELISATION START]");
	} else {

		MPI_Status status_recv;
		Branch branch_recv;
    
    	// recv work
    	MPI_Recv(&branch_recv, sizeof(Branch), MPI_BYTE, 0, TAG_WORK, MPI_COMM_WORLD, &status_recv);
		std::atomic<unsigned short> best_ub = branch_recv.ub;

		queue.push(std::move(branch_recv));
		
		int best_ub_unchanged_count = 0; //do we need this?
	}

	// OpenMP Parallel Region
	/*
	idea is assign specific threads to specific tasks, in particular the first three threads are used for 
	gathering best_ub, checking for termination and listening for requests, than the other threads are used for
		-one to keep popping from the queue, work_stealing, generate omp_tasks(method create_task) and add new branches to the queue
		-others to compute the omp_tasks
	*/
	#pragma omp parallel shared(best_ub, queue, queue_mutex, active_tasks, max_tasks)
	{
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
		// Only worker processes go in here.
		if (my_rank != 0) {
			//worker thread used for listening for requests about work stealing
			if (tid == 2) {
				thread_2_listen_for_requests(queue_mutex, queue);
			}

			Branch current;
			std::unique_ptr<Graph> current_G;
			int current_lb;
			unsigned short current_ub;
			unsigned int u,v;
			std::vector<Branch> new_branches(2);

			#pragma omp single nowait shared(queue, queue_mutex, best_ub, active_tasks, max_tasks)
			{
				MPI_Status status_work_steal;

				//TODO:change true with flag terminate (is false until we receive a signal from terminator)
				while (true) { // keep dequeuing until queue is empty or we dont have other threads available(or too much tasks generated)
					bool has_work = false;
					
					{
						std::lock_guard<std::mutex> lock(queue_mutex);
						if (!queue.empty()) {
							current = std::move(const_cast<Branch&>(queue.top()));
							queue.pop();
							has_work = true;
						}
					}

					//TODO: if response is 0 ask to others(keep asking to others? limit the number of requests?)
					// Work Stealing
					if (!has_work) {

						int target_worker = (my_rank + rand() % (p - 1) + 1) % p; //check if it is correct

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

					Log_par("Processing node: lb = " + std::to_string(current_lb) +
						", ub = " + std::to_string(current_ub),
		    			current.depth);

					// find solution
					if (current_lb == current_ub) {
						MPI_Send(&current_ub, 1, MPI_INT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD); //check if it is correct
						Log_par("[FOUND] Chromatic number found: " +
							std::to_string(current_lb),
							current.depth);
						Log_par("========== END ==========", 0);
						continue;
					}

					// Prune
					if (current_lb >= best_ub.load()){
						Log_par("[PRUNE] Branch pruned at depth " +
							std::to_string(current.depth) +
							": lb = " + std::to_string(current_lb) +
							" >= best_ub = " + std::to_string(best_ub),
							current.depth);
						continue;
					}

					auto type = _branching_strat.PairType::DontCare;
					std::pair<unsigned int, unsigned int> vertices = _branching_strat.ChooseVertices(*current_G, type);
					u = vertices.first;
					v = vertices.second;
					Log_par("Branching on vertices: u = " + std::to_string(u) +
							", v = " + std::to_string(v),
		    				current.depth, true);

					if (u == -1 || v == -1) {
						int chromatic_number = current_G->GetNumVertices();
						MPI_Send(&chromatic_number, 1, MPI_INT, 0, TAG_SOLUTION_FOUND, MPI_COMM_WORLD); //check if it is correct
						Log_par("Graph is complete. Chromatic number = " +
							std::to_string(chromatic_number),
							current.depth);
						Log_par("========== END ==========", 0);
						continue;
					}

					//TODO: check for other approaches to avoid generating too much tasks(this will wait untill all child tasks are completed)
					//avoid generating too much tasks
					while (active_tasks > max_tasks) {
						#pragma omp taskwait
					}

					//generate tasks and update queue
					#pragma omp task firstprivate(current_G, u, v)
					{
						create_task(active_tasks, current_G.get(), u, v, _clique_strat, _color_strat, new_branches, 1, best_ub, current.depth);
						create_task(active_tasks, current_G.get(), u, v, _clique_strat, _color_strat, new_branches, 2, best_ub, current.depth);

						{
							std::lock_guard<std::mutex> lock(queue_mutex);
							for (auto& branch : new_branches) {
								if (branch.g) queue.push(std::move(branch));
							}
						}
					}

					// Update local sbest_ub
					unsigned short previous_best_ub = best_ub.load();
					best_ub.store(std::min({previous_best_ub, new_branches[0].ub, new_branches[1].ub}));
					Log_par("[UPDATE] Updated best_ub: " + std::to_string(best_ub),
		    		current.depth);
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

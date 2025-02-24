# Graph Color Number

## Overview
This project implements a high-performance algorithm to find the chromatic number of a graph. The algorithm uses a parallel branch & bound approach with heuristics to determine lower and upper bounds for the chromatic number at each branch.

This project is a solution to the student challenge for EuroHPC Summit 2025.

## Prerequisites
- MPI (Message Passing Interface) library
- OpenMP (Open Multi-Processing) library
- C++17 compatible compiler

## Installation
1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/graph-color-number.git
    cd graph-color-number
    ```

2. Build the project:
    ```sh
    mkdir build
    cd build
    cmake ..
    make
    ```

## Usage
To run the script, use the following command:
```sh
mpirun -np <number_of_processes> ./src/scripts/run_instance <file_name> [timeout]
```
- `<number_of_processes>`: Number of MPI processes to use.
- `<file_name>`: Name of the graph file located in the `graphs_instances` directory.
- `[timeout]`: (Optional) Timeout in seconds. Default is 60 seconds.

Example:
```sh
mpirun -np 4 ./src/scripts/run_instance anna.col 120
```

**Note:** It is recommended to use OpenMPI/4.1.4-GCC-11.3.0 and CMake/3.23.1-GCCcore-11.3.0.

## Running on a Supercomputer (Vega) with Slurm
### Load Necessary Modules
```sh
module load OpenMPI/4.1.4-GCC-11.3.0
module load CMake/3.23.1-GCCcore-11.3.0
```

### Example Slurm Script
Create a file named `run_instance_job.slurm` with the following content:
```bash
#!/bin/bash
#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name
#SBATCH --nodes=64                             # Number of nodes
#SBATCH --ntasks-per-node=8                    # MPI tasks per node
#SBATCH --cpus-per-task=16                     # Cores per MPI task
#SBATCH --time=01:00:00                        # Time limit (HH:MM:SS)
#SBATCH --partition=cpu                        
#SBATCH --output=output_mpi.txt                # Standard output file
#SBATCH --error=error_mpi.txt                  # Standard error file

# Run the MPI program
cd build/src/scripts/
srun run_instance le450_15a.col 1800
```

Submit the job using:
```sh
sbatch run_job.slurm
```

## Expected Results
The script compares the computed chromatic number with the expected results stored in `expected_chi.txt`. If the computed result does not match the expected result, an error message is displayed.

## Logs
Logs are generated for each MPI process and stored in the `logs` directory. The log files are named `log_<rank>.txt`, where `<rank>` is the MPI process rank. It contains detailed information on each branch's intermediate results (lower and upper bounds). 

## run_instance.cpp Script Details
- The script reads a graph file and initializes the MPI environment.
- It loads the expected chromatic number from `expected_chi.txt`.
- The graph is processed using the parallel branch & bound algorithm to find the chromatic number.
- The computed chromatic number is compared with the expected result.
- The script logs the execution details and results.

## Troubleshooting
- Ensure that the graph files are placed in the `graphs_instances` directory.
- Verify that the `expected_chi.txt` file contains the expected chromatic numbers for the graph files.
- Check the log files in the `logs` directory for detailed execution information.



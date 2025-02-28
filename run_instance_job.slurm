#!/bin/bash
#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name
#SBATCH --nodes=32                            # Number of nodes
#SBATCH --ntasks-per-node=8                    # MPI tasks per node
#SBATCH --cpus-per-task=4                      # Cores per MPI task (increase if CPU-BOUND)
#SBATCH --time=00:10:00                        # Time limit (HH:MM:SS)
#SBATCH --partition=cpu                        # Replace with actual partition
#SBATCH --output=output_mpi.txt
#SBATCH --error=error_mpi.txt

# Extract the base name of the graph file (without extension)
GRAPH_FILE="le450_15a.col"

# Use the base name for output and error files

# Run the MPI program
cd build/src/scripts/
srun run_instance "$GRAPH_FILE" --timeout=500 --sol_gather_period=40 --balanced=1 --color_strategy=2 --output="${GRAPH_FILE}_output.txt" --logging=0

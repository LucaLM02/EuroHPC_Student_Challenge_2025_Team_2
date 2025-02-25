#!/bin/bash
#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name
#SBATCH --nodes=64                        # Number of nodes
#SBATCH --ntasks-per-node=8              # MPI tasks per node
#SBATCH --cpus-per-task=16               # Cores per MPI task
#SBATCH --time=03:00:00                  # Time limit (HH:MM:SS)
#SBATCH --partition=cpu                  # Replace with actual partition
#SBATCH --output=output_mpi_%A_%a.txt    # Output file with job array ID
#SBATCH --error=error_mpi_%A_%a.txt      # Error file with job array ID

# Array of graph files in the directory
GRAPH_FILES=(build/src/scripts/graph_instances/*.col)

# Get the file for the current job array index
GRAPH_FILE="${GRAPH_FILES[$SLURM_ARRAY_TASK_ID]}"

# Run the MPI program
cd build/src/scripts/
srun run_instance "$GRAPH_FILE" 100000 20 


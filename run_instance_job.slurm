#!/bin/bash
#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name
#SBATCH --nodes=64                        # Number of nodes
#SBATCH --ntasks-per-node=8              # MPI tasks per node
#SBATCH --cpus-per-task=16               # Cores per MPI task (increase if CPU-BOUND)
#SBATCH --time=03:00:00                  # Time limit (HH:MM:SS)
#SBATCH --partition=cpu                  # Replace with actual partition
#SBATCH --output=output_mpi.txt           # Standard output file
#SBATCH --error=error_mpi.txt             # Standard error file

# Run the MPI program
cd build/src/scripts/
srun run_instance le450_15a.col 100000 


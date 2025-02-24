#!/bin/bash
#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name
#SBATCH --nodes=8                        # Number of nodes
#SBATCH --ntasks-per-node=1              # MPI tasks per node
#SBATCH --cpus-per-task=5               # Cores per MPI task
#SBATCH --time=00:10:00                  # Time limit (HH:MM:SS)
#SBATCH --partition=cpu                  # Replace with actual partition
#SBATCH --output=output_mpi.txt           # Standard output file
#SBATCH --error=error_mpi.txt             # Standard error file

# Check network interfaces and InfiniBand status
echo "Checking available network interfaces:"
ip a

echo "Checking InfiniBand status:"
ibstat

# Run the MPI program
cd build/tests/branch_n_bound_par/
srun test_b_n_b_par le450_15a.col 300 1


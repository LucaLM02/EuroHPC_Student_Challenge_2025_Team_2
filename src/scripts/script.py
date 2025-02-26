import time
import os
import subprocess

# Paths
SLURM_FILE = "run.slurm"
OUTPUT_FILE = "chromatic_number.txt"
TIME_FILE = "time.txt"

# Define chromatic numbers (example mapping)
chromatic_numbers = {
        "anna.col":         [11,    1],
        "fpsol2.i.1.col":   [65,    1],
        "fpsol2.i.3.col":   [30,    1],

        "inithx.i.2.col":   [31,    2],
        "le450_5c.col":     [5,     2],
        "queen6_6.col":     [7,     4],
        "queen9_9.col":     [10,    4],

        "le450_5a.col":     [5,     8],
        "queen13_13.col":   [13,    8],
        "le450_15d.col":    [15,    8]
}  # Update with real values

# Execution time (adjust as needed)
EXECUTION_TIME = 10000

# Iterate over all files in the graphs directory
for filename in chromatic_numbers.keys():
    file_path = filename

    # Skip if the file has no associated chromatic number

    chromatic_number = chromatic_numbers[filename][0]

    for color_strategy in range(2):
        for isBalanced in range(2):

            output_file = filename.split(".")[0] + "_output_" + str(isBalanced) + "_" + str(color_strategy) + ".txt" 

            # Generate the SLURM script
            slurm_script = f"""#!/bin/bash
        #SBATCH --job-name=team2_mpi_branch_n_bound
        #SBATCH --nodes={chromatic_numbers[filename][1]}
        #SBATCH --ntasks-per-node=8
        #SBATCH --cpus-per-task=16
        #SBATCH --time=03:00:00
        #SBATCH --partition=cpu
        #SBATCH --output=output_mpi.txt
        #SBATCH --error=error_mpi.txt

        srun test_graph {file_path} {chromatic_number} {EXECUTION_TIME} {output_file} {color_strategy} {isBalanced}
        """

            # Write the SLURM script
            with open(SLURM_FILE, "w") as f:
                f.write(slurm_script)

            print(f"Submitting job for {filename}...")

            # Submit the SLURM job and capture job ID
            submit_output = subprocess.run(["sbatch", SLURM_FILE], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

            
            job_id = None
            for word in submit_output.stdout.split():
                if word.isdigit():
                    job_id = word
                    break
            
            if not job_id:
                print(f"Failed to submit job for {filename}")
                continue

            print(f"Job {job_id} submitted for {filename}, waiting for completion...")

            # Efficiently wait for job completion
            while True:
                check_job = subprocess.run(["squeue", "-j", job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                if job_id not in check_job.stdout:
                    break  # Job is no longer in queue
                time.sleep(10)  # Wait 10 seconds before checking again

            print(f"Job {job_id} completed. Reading output...")




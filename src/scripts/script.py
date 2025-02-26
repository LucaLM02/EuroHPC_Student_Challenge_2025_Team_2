import time
import os
import subprocess

# Paths
SLURM_FILE = "run.slurm"
OUTPUT_FILE = "chromatic_number.txt"
TIME_FILE = "time.txt"

# Define chromatic numbers (example mapping)
chromatic_numbers = {
        "anna.col": 11,
        "david.col": 11,
        "fpsol2.i.1.col": 65,
        "fpsol2.i.2.col": 30,
        "fpsol2.i.3.col": 30,
        "games120.col": 9,
        "homer.col": 13,
        "huck.col": 11,
        "inithx.i.1.col": 54,
        "inithx.i.2.col": 31,
        "inithx.i.3.col": 31,
        "jean.col": 10,
        "le450_5a.col": 5,
        "le450_5b.col": 5,
        "le450_5c.col": 5,
        "le450_5d.col": 5,
        "le450_15a.col": 15,
        "le450_15b.col": 15,
        "le450_15c.col": 15,
        "le450_15d.col": 15,
        "le450_25a.col": 25,
        "le450_25b.col": 25,
        "le450_25c.col": 25,
        "le450_25d.col": 25,
        "miles250.col": 8,
        "miles500.col": 20,
        "miles750.col": 31,
        "miles1000.col": 42,
        "miles1500.col": 73,
        "mulsol.i.1.col": 49,
        "mulsol.i.2.col": 31,
        "mulsol.i.3.col": 31,
        "mulsol.i.4.col": 31,
        "mulsol.i.5.col": 31,
        "myciel3.col": 4,
        "myciel4.col": 5,
        "myciel5.col": 6,
        "myciel6.col": 7,
        "myciel7.col": 8,
        "queen5_5.col": 5,
        "queen6_6.col": 7,
        "queen7_7.col": 7,
        "queen8_8.col": 9,
        "queen8_12.col": 12,
        "queen9_9.col": 10,
        "queen11_11.col": 11,
        "queen13_13.col": 13,
        "zeroin.i.1.col": 49,
        "zeroin.i.2.col": 30,
        "zeroin.i.3.col": 30
}  # Update with real values

# Execution time (adjust as needed)
EXECUTION_TIME = 180 

# Iterate over all files in the graphs directory
for filename in chromatic_numbers.keys():
    file_path = filename

    # Skip if the file has no associated chromatic number

    chromatic_number = chromatic_numbers[filename]

    # Generate the SLURM script
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=team2_mpi_branch_n_bound
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --time=03:00:00
#SBATCH --partition=cpu
#SBATCH --output=output_mpi.txt
#SBATCH --error=error_mpi.txt

srun test_b_b_n_b_par {file_path} {chromatic_number} {EXECUTION_TIME}
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

    # Read the result from output.txt
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, "r") as f:
            result = f.read().strip()
        if os.path.exists(TIME_FILE):
            with open(TIME_FILE, "r") as ff:
                time_result = ff.read().strip()
            print(f"Chromatic number of {filename} is {result} (expected {chromatic_number}) in {time_result}")
        else:
            print(f"Warning: {TIME_FILE} not found for {filename}")
    else:
        print(f"Warning: {OUTPUT_FILE} not found for {filename}")


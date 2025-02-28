#include <string>
#include <unordered_map>
#include <iostream>
#include <cstdlib> // For std::stoi
#include <fstream>
#include <filesystem>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

int main(int argc, char* argv[]) {
    std::unordered_map<std::string, int> expectedEasyResults;
    std::unordered_map<std::string, int> expectedMediumResults;
    std::unordered_map<std::string, int> expectedHardResults;
    expectedEasyResults.insert({"anna.col", 11});
    expectedEasyResults.insert({"david.col", 11});
    expectedEasyResults.insert({"fpsol2.i.1.col", 65});
    expectedEasyResults.insert({"fpsol2.i.2.col", 30});
    expectedEasyResults.insert({"fpsol2.i.3.col", 30});
    expectedEasyResults.insert({"games120.col", 9});
    expectedEasyResults.insert({"homer.col", 13});
    expectedEasyResults.insert({"huck.col", 11});
    expectedEasyResults.insert({"inithx.i.1.col", 54});
    expectedEasyResults.insert({"inithx.i.2.col", 31});
    expectedEasyResults.insert({"inithx.i.3.col", 31});
    expectedEasyResults.insert({"jean.col", 10});
    expectedEasyResults.insert({"miles250.col", 8});
    expectedEasyResults.insert({"miles500.col", 20});
    expectedEasyResults.insert({"miles750.col", 31});
    expectedEasyResults.insert({"miles1000.col", 42});
    expectedEasyResults.insert({"miles1500.col", 73});
    expectedEasyResults.insert({"myciel3.col", 4});
    expectedEasyResults.insert({"myciel4.col", 5});
    expectedEasyResults.insert({"myciel5.col", 6});
    expectedEasyResults.insert({"myciel6.col", 7});
    expectedEasyResults.insert({"myciel7.col", 8});
    expectedEasyResults.insert({"zeroin.i.1.col", 49});
    expectedEasyResults.insert({"zeroin.i.2.col", 30});
    expectedEasyResults.insert({"zeroin.i.3.col", 30});

    expectedMediumResults.insert({"queen5_5.col", 5});
    expectedMediumResults.insert({"queen6_6.col", 7});
    expectedMediumResults.insert({"queen7_7.col", 7});

    expectedHardResults.insert({"queen8_8.col", 9});
    expectedHardResults.insert({"queen8_12.col", 12});
    expectedHardResults.insert({"queen9_9.col", 10});
    expectedHardResults.insert({"queen11_11.col", 11});
    expectedHardResults.insert({"queen13_13.col", 13});
    expectedHardResults.insert({"le450_5a.col", 5});
    expectedHardResults.insert({"le450_5b.col", 5});
    expectedHardResults.insert({"le450_5c.col", 5});
    expectedHardResults.insert({"le450_5d.col", 5});
    expectedHardResults.insert({"le450_15a.col", 15});
    expectedHardResults.insert({"le450_15b.col", 15});
    expectedHardResults.insert({"le450_15c.col", 15});
    expectedHardResults.insert({"le450_15d.col", 15});
    expectedHardResults.insert({"le450_25a.col", 25});
    expectedHardResults.insert({"le450_25b.col", 25});
    expectedHardResults.insert({"le450_25c.col", 25});
    expectedHardResults.insert({"le450_25d.col", 25});
    expectedHardResults.insert({"mulsol.i.1.col", 49});
    expectedHardResults.insert({"mulsol.i.3.col", 31});
    expectedHardResults.insert({"mulsol.i.4.col", 31});
    expectedHardResults.insert({"mulsol.i.5.col", 31});

    bool run_easy = false;
    bool run_medium = false;
    bool run_hard = false;

    // Parse optional arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--run_easy") {
            run_easy = true;
        } else if (arg == "--run_medium") {
            run_medium = true;
        } else if (arg== "--run_hard") {
            run_hard = true;
        } else {
            std::cerr << "Error: Invalid argument format " << arg << ".\n";
            return 1;
        }
    }

    if ( run_easy ) {
        for ( auto file_chi_pair : expectedEasyResults ) {
            std::string command;
            std::string output_file = "run_all_instances.slurm";
            std::ofstream ofs(output_file);
            
            ofs << "#!/bin/bash " << std::endl
                << "#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name" << std::endl 
                << "#SBATCH --nodes=1                        # Number of nodes" << std::endl
                << "#SBATCH --ntasks-per-node=8              # MPI tasks per node" << std::endl
                << "#SBATCH --cpus-per-task=4               # Cores per MPI task (increase if CPU-BOUND)" << std::endl
                << "#SBATCH --time=00:15:00                  # Time limit (HH:MM:SS)" << std::endl
                << "#SBATCH --partition=cpu                  # Replace with actual partition" << std::endl
                << "#SBATCH --output=output_mpi.txt           # Standard output file" << std::endl
                << "#SBATCH --error=error_mpi.txt             # Standard error file" << std::endl << std::endl;
            
            std::string output = file_chi_pair.first;
            output = output.substr(0, output.size() - 3) + "_output.txt";
            ofs << "srun run_instance " << file_chi_pair.first << " --timeout=600 --sol_gather_period=20 --balanced=1 --output=" 
                << output << std::endl;

            command = "sbatch " + output_file;
            std::cout << exec(command.c_str()) << std::endl;
        }
    }

    if ( run_medium ) {
        for ( auto file_chi_pair : expectedMediumResults ) {
            std::string command;
            std::string output_file = "run_all_instances.slurm";
            std::ofstream ofs(output_file);
            
            ofs << "#!/bin/bash " << std::endl
                << "#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name" << std::endl 
                << "#SBATCH --nodes=8                        # Number of nodes" << std::endl
                << "#SBATCH --ntasks-per-node=8              # MPI tasks per node" << std::endl
                << "#SBATCH --cpus-per-task=4               # Cores per MPI task (increase if CPU-BOUND)" << std::endl
                << "#SBATCH --time=04:00:00                  # Time limit (HH:MM:SS)" << std::endl
                << "#SBATCH --partition=cpu                  # Replace with actual partition" << std::endl
                << "#SBATCH --output=output_mpi.txt           # Standard output file" << std::endl
                << "#SBATCH --error=error_mpi.txt             # Standard error file" << std::endl << std::endl;
            
            std::string output = file_chi_pair.first;
            output = output.substr(0, output.size() - 3) + "_output.txt";
            ofs << "srun run_instance " << file_chi_pair.first << " --timeout=10000 --sol_gather_period=20 --balanced=1 --output=" 
                << output << std::endl;


            command = "sbatch " + output_file;
            std::cout << exec(command.c_str()) << std::endl;
        }
    }

    if ( run_hard ) {
        for ( auto file_chi_pair : expectedHardResults ) {
            std::string command;
            std::string output_file = "run_all_instances.slurm";
            std::ofstream ofs(output_file);
            
            ofs << "#!/bin/bash " << std::endl
                << "#SBATCH --job-name=team2_mpi_branch_n_bound    # Job name" << std::endl 
                << "#SBATCH --nodes=64                        # Number of nodes" << std::endl
                << "#SBATCH --ntasks-per-node=8              # MPI tasks per node" << std::endl
                << "#SBATCH --cpus-per-task=4               # Cores per MPI task (increase if CPU-BOUND)" << std::endl
                << "#SBATCH --time=04:00:00                  # Time limit (HH:MM:SS)" << std::endl
                << "#SBATCH --partition=cpu                  # Replace with actual partition" << std::endl
                << "#SBATCH --output=output_mpi.txt           # Standard output file" << std::endl
                << "#SBATCH --error=error_mpi.txt             # Standard error file" << std::endl << std::endl;
            
            std::string output = file_chi_pair.first;
            output = output.substr(0, output.size() - 3) + "_output.txt";
            ofs << "srun run_instance " << file_chi_pair.first << " --timeout=10000 --sol_gather_period=20 --balanced=1 --output=" 
                << output << std::endl;

            command = "sbatch " + output_file;
            std::cout << exec(command.c_str()) << std::endl;
        }
    }

}

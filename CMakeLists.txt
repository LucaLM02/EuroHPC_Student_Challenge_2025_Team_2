cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Enable testing (useful for CTest)
# enable_testing()

# Include subdirectories
<<<<<<< Updated upstream
add_subdirectory(src)                # Build source files first
add_subdirectory(tests/dimacs_graph) # Build constraints test
add_subdirectory(tests/common)       # Build test utilities
add_subdirectory(tests/color)        # Build color test
=======
add_subdirectory(src)                   # Build source files first
add_subdirectory(tests/common)          # Build test utilities
add_subdirectory(tests/color)           # Build color test
add_subdirectory(tests/dimacs_graph)    # Build dimacs graph test
add_subdirectory(tests/csr_graph)       # Build csr graph test
add_subdirectory(tests/branch_n_bound)  # Build branch_n_bound test
>>>>>>> Stashed changes

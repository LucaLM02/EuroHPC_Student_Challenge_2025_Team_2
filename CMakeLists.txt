cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Enable testing (useful for CTest)
# enable_testing()

# Include subdirectories
add_subdirectory(src)                   # Build source files first
add_subdirectory(tests/clique_strategy) # Build dsatur color test
add_subdirectory(tests/common)          # Build test utilities
add_subdirectory(tests/color)           # Build color test
add_subdirectory(tests/dsatur_color)    # Build dsatur color test
add_subdirectory(tests/dimacs_graph)    # Build dimacs graph test
add_subdirectory(tests/csr_graph)       # Build csr graphc test
add_subdirectory(tests/branch_n_bound)  # Build branch_n_bound test

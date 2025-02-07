cmake_minimum_required(VERSION 3.10)
project(test_clique_strategy)

set(CMAKE_CXX_STANDARD 20)

# Add OpenMP support
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Explicitly tell CMake to include the "include" directory globally
include_directories(${PROJECT_SOURCE_DIR}/src)

# Ensure headers are accessible in all packages (src, include, tests)
include_directories(${PROJECT_SOURCE_DIR}/src/clique)
include_directories(${PROJECT_SOURCE_DIR}/src/base)
include_directories(${PROJECT_SOURCE_DIR}/tests)

# Add the directory where fastwclq.hpp is located
include_directories(${PROJECT_SOURCE_DIR}/src/clique)

# Add executable including main.cpp
add_executable(test_clique_strategy
        src/clique/fastwclq.cpp
        src/base/graph_clique.cpp
        tests/clique_strategy/test.cpp)

# Ensure target includes all directories for global header access
target_include_directories(test_clique_strategy PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(test_clique_strategy PRIVATE ${PROJECT_SOURCE_DIR}/tests)
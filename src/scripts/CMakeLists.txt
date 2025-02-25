find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER mpicxx)
include_directories(${MPI_INCLUDE_DIR})

SET(GCC_MY_COMPILE_FLAGS "-g -std=c++20") #"-g3 -std=c++20")
SET(GCC_MY_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_MY_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_MY_LINK_FLAGS}")

add_executable(run_instance run_instance.cpp)

# Link test_color executable with the main library and common test utilities
target_link_libraries(run_instance PRIVATE chromatic_number test_common)

# Include necessary headers
target_include_directories(run_instance PRIVATE 
    ${CMAKE_SOURCE_DIR}/src)

# Copy necessary files and directories to the build/src/scripts directory
file(COPY ${CMAKE_SOURCE_DIR}/src/scripts/expected_chi.txt DESTINATION ${CMAKE_BINARY_DIR}/src/scripts)
file(COPY ${CMAKE_SOURCE_DIR}/src/scripts/graphs_instances DESTINATION ${CMAKE_BINARY_DIR}/src/scripts)
file(COPY ${CMAKE_SOURCE_DIR}/src/scripts/logs DESTINATION ${CMAKE_BINARY_DIR}/src/scripts)

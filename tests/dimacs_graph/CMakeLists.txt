SET(GCC_MY_COMPILE_FLAGS "-g -std=c++20")  #"-g3 -std=c++20")
SET(GCC_MY_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_MY_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_MY_LINK_FLAGS}")

add_executable(test_dimacs test.cpp)

# Link test_color executable with the main library and common test utilities
target_link_libraries(test_dimacs PRIVATE chromatic_number test_common)

# Include necessary headers
target_include_directories(test_dimacs PRIVATE 
    ${CMAKE_SOURCE_DIR}/src 
    ${CMAKE_SOURCE_DIR}/tests/common)
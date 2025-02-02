SET(GCC_MY_COMPILE_FLAGS "-g -std=c++20")  #"-g3 -std=c++20")
SET(GCC_MY_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_MY_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_MY_LINK_FLAGS}")

add_library(test_common STATIC test_common.cpp)

# Make headers available to test executables
target_include_directories(test_common PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}  # Includes tests/common/
    ${CMAKE_SOURCE_DIR}/src      # Includes src/ 
    ${CMAKE_SOURCE_DIR}/src/base # Includes src/base
)

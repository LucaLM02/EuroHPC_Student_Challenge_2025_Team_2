SET(GCC_MY_COMPILE_FLAGS "-g -std=c++20")  #"-g3 -std=c++20")
SET(GCC_MY_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_MY_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_MY_LINK_FLAGS}")

# Find all source files in src/ and src/base/
file(GLOB SRC_FILES common.cpp *.cpp color/*.cpp base/*.cpp branching/*.cpp branch_n_bound/branch_n_bound_seq.cpp clique/*.cpp)

# Create a static library from all source files
add_library(chromatic_number STATIC ${SRC_FILES})

# Add both src/ and src/base/ to the include directories
target_include_directories(chromatic_number PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}                  # Includes src/
    ${CMAKE_CURRENT_SOURCE_DIR}/base             # Includes src/base/
    ${CMAKE_CURRENT_SOURCE_DIR}/color            # Includes src/color/
    ${CMAKE_CURRENT_SOURCE_DIR}/branching        # Includes src/branching/
    ${CMAKE_CURRENT_SOURCE_DIR}/branch_n_bound   # Includes src/branch_n_bound/
    ${CMAKE_CURRENT_SOURCE_DIR}/clique           # Includes src/clique/
)
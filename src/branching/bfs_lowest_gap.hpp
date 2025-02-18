#ifndef BFS_LOWEST_GAP_HPP
#define BFS_LOWEST_GAP_HPP

#include <vector>
#include <queue>
#include <set>
#include "csr_graph.hpp"
#include "branching_strategy.hpp"
#include "graph.hpp"
#include "fastwclq.hpp"



// Function prototypes
int selectBranchingVertex(const CSRGraph &graph, const std::vector<int> &coloring);

#endif  // BFS_LOWEST_GAP_HPP
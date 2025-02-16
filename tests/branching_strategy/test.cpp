#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include "branching_strategy.hpp"
#include "csr_graph.hpp"

#include "test_common.hpp"



int main() {
    const std::string file_name = "10_vertices_graph.col";

    Graph& graph = *CSRGraph::LoadFromDimacs(file_name);

    NeighboursBranchingStrategy branching_strategy;

    auto [v, u] = branching_strategy.ChooseVertices(graph);

    std::cout << "Vertices with the greatest number of common neighbours are " 
              << v << " and " << u << std::endl;

    return 0;
}
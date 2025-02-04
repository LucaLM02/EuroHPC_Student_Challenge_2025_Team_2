#include <iostream>
#include <memory>

#include "graph.hpp"
#include "dimacs_graph.hpp"
#include "csr_graph.hpp"
#include "color.hpp"
#include "test_common.hpp"

void test_graph(Graph& graph) {
    
    GreedyColorStrategy color_strategy;

    std::vector<unsigned short> coloring;
    unsigned short max_k;

    color_strategy.Color(graph, coloring, max_k);
    
    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(coloring) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
}

void test_dimacs_graph(const std::string& file_name) {

    DimacsGraph* graph = DimacsGraph::LoadFromDimacs(file_name);
    test_graph(*graph);
}

void test_csr_graph(const std::string& file_name) {
    CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);
    test_graph(*graph);
}

int main() {
    const std::string file_name = "10_vertices_graph.clq";

    std::cout << "-- COLORING DIMACS GRAPH --" << std::endl;
    test_dimacs_graph(file_name);
    std::cout << std::endl;

    std::cout << "-- COLORING DIMACS GRAPH --" << std::endl;
    test_csr_graph(file_name);

    return 0;
}
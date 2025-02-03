#include <iostream>
#include <memory>

#include "graph.hpp"
#include "dimacs_graph.hpp"
#include "color.hpp"
#include "test_common.hpp"

void test_graph(Graph&& graph) {
    
    GreedyColorStrategy color_strategy;

    std::vector<unsigned short> coloring;
    unsigned short max_k;

    color_strategy.Color(graph, coloring, max_k);
    
    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(coloring) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
}


int main() {
    std::string file_name = "10_vertices_graph";
    Dimacs dimacs;
    dimacs.load(file_name.c_str());

    test_graph(DimacsGraph(dimacs));
    return 0;
}
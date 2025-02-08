#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include "graph.hpp"
#include "dimacs_graph.hpp"
#include "csr_graph.hpp"

#include "color.hpp"

#include "test_common.hpp"

void test_graph(Graph& graph) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        GreedyColorStrategy color_strategy;

        unsigned short max_k;

        color_strategy.Color(graph, max_k);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    std::cout << "Time to color a graph with " << graph.GetNumVertices() << "vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(graph.GetColoring()) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
    std::cout << "Is valid: " << TestFunctions::CheckColoring(graph) << std::endl;

}

void test_dimacs_graph(const std::string& file_name) {

    DimacsGraph* graph = DimacsGraph::LoadFromDimacs(file_name);
    test_graph(*graph);
}

void test_csr_graph(const std::string& file_name) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    CSRGraph& graph = *CSRGraph::LoadFromDimacs(file_name);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    std::cout << "Time to load a graph with " << graph.GetNumVertices() << "vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    graph.SortByDegree(false);

    test_graph(graph);

    std::cout << " ======================== COLORING WITH FEEDBACK ORDER ========================" << std::endl;

    test_graph(graph);
}


int main() {
    const std::string file_name = "queen15_15.col";

    std::cout << "-- COLORING DIMACS GRAPH --" << std::endl;
    //test_dimacs_graph(file_name);
    std::cout << std::endl;

    std::cout << "-- COLORING CSR GRAPH --" << std::endl;
    test_csr_graph(file_name);

    return 0;
}
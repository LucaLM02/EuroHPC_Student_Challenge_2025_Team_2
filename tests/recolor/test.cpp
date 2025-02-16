#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include "graph.hpp"
#include "csr_graph.hpp"

#include "dsatur_color.hpp"
#include "recolor.hpp"

#include "test_common.hpp"

void test_graph(Graph& graph) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        DSaturColorStrategy color_strategy;

        unsigned short max_k;

        color_strategy.Color(graph, max_k);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    std::cout << "Time to color a graph with " << graph.GetNumVertices() << " vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(graph.GetColoring()) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
    std::cout << "Is valid: " << TestFunctions::CheckColoring(graph) << std::endl;

    {
    begin = std::chrono::steady_clock::now();

        GreedySwapRecolorStrategy recolor_strategy;

        if ( recolor_strategy.Recolor(graph) ) {
            std::cout << "Successfully recolored the graph" << std::endl;
        } else {
            std::cout << "Unable to recolor the graph" << std::endl;
        }

    end = std::chrono::steady_clock::now();

    std::cout << "Time to recolor a graph with " << graph.GetNumVertices() << " vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    std::vector<unsigned short> coloring = graph.GetColoring();
    std::sort(coloring.begin(), coloring.end());
    max_k = coloring[coloring.size()-1];

    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(graph.GetColoring()) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
    std::cout << "Is valid: " << TestFunctions::CheckColoring(graph) << std::endl;
    }

    {
        begin = std::chrono::steady_clock::now();

            graph.SortByColor(false);
            GreedyColorStrategy color_strategy;

            color_strategy.Color(graph, max_k);

        end = std::chrono::steady_clock::now();

        std::cout << "Time to color a graph with " << graph.GetNumVertices() << " vertices and " 
                << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;


        std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
        std::cout << "Coloring: " << TestFunctions::VecToString(graph.GetColoring()) << std::endl;
        std::cout << "Max k:    " << max_k << std::endl;
        std::cout << "Is valid: " << TestFunctions::CheckColoring(graph) << std::endl;
    }

    {
    begin = std::chrono::steady_clock::now();

        GreedySwapRecolorStrategy recolor_strategy;

        if ( recolor_strategy.Recolor(graph) ) {
            std::cout << "Successfully recolored the graph" << std::endl;
        } else {
            std::cout << "Unable to recolor the graph" << std::endl;
        }

    end = std::chrono::steady_clock::now();

    std::cout << "Time to recolor a graph with " << graph.GetNumVertices() << " vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    std::vector<unsigned short> coloring = graph.GetColoring();
    std::sort(coloring.begin(), coloring.end());
    max_k = coloring[coloring.size()-1];

    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(graph.GetColoring()) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
    std::cout << "Is valid: " << TestFunctions::CheckColoring(graph) << std::endl;
    }

}

void test_csr_graph(const std::string& file_name) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        CSRGraph& graph = *CSRGraph::LoadFromDimacs(file_name);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    std::cout << "Time to load a graph with " << graph.GetNumVertices() << "vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    test_graph(graph);

}


int main() {
    const std::string file_name = "school1.col";

    std::cout << "-- COLORING DIMACS GRAPH --" << std::endl;
    //test_dimacs_graph(file_name);
    std::cout << std::endl;

    std::cout << "-- COLORING CSR GRAPH --" << std::endl;
    test_csr_graph(file_name);

    return 0;
}
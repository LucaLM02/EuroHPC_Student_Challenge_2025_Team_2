#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include "graph.hpp"
#include "csr_graph.hpp"

#include "color.hpp"
#include "dsatur_color.hpp"

#include "test_common.hpp"

void basic_dsatur_list_test(const DSaturList& list, int indentation=1) {
    std::string indentation_string;

    for ( int i=0; i<indentation; i++) {
        indentation_string.append("  ");
    }


    int highest_sat_degree = list.GetHighestSatDegree();

    std::cout << indentation_string << "Highest vertex: " << list.GetHighestVertex() 
              << " with sat degree: " << highest_sat_degree << std::endl;
    std::cout << indentation_string << "Lowest vertex: " << list.GetLowestVertex() 
              << " with sat degree: " << list.GetLowestSatDegree() << std::endl;
    std::cout << indentation_string << "Is empty: " << (list.IsEmpty() ? "true" : "false") 
              << std::endl;
    for (int i = 0; i <= highest_sat_degree; i++) {
        const DSaturItem* item = list[i];
        std::cout << indentation_string << "sat=" << i;
        while ( item != nullptr ) {
            std::cout << " (" << item->vertex << "," << item->degree <<")";
            item = item->next;
        }
        std::cout << std::endl;
    }

}

void test_dsatur_list(const std::string& file_name) {
    CSRGraph& graph = *CSRGraph::LoadFromDimacs(file_name);

    DSaturList list(graph);
    std::vector<int> neighbours;

    std::cout << " ==================== DSaturList Creation ==================== " << std::endl;
    basic_dsatur_list_test(list);

    int max_k = 0;
    int selected_vertex;
    unsigned short selected_color;
    std::vector<unsigned short> coloring(graph.GetHighestVertex()+1);

    std::cout << " ======================= Graph Coloring ====================== " << std::endl;
    while ( !list.IsEmpty() ) {
        // retrieves the element with highest degree among the ones with highest 
        // saturation degree
        selected_vertex = list.PopHighestVertex();

        // colors the vertex in a greedy fashion
        selected_color  = GreedyFindColor(graph, selected_vertex, coloring, max_k + 1);

        std::cout << "   Selected vertex: " << selected_vertex 
                  << " Selected color: "  << selected_color << std::endl;
        coloring[selected_vertex] = selected_color;
        std::cout << "   Right after coloring: " << std::endl;
        basic_dsatur_list_test(list, 2);

        // updates the maximum color
        if ( selected_color > max_k ) {
            max_k = selected_color;
        }

        // updating the saturation degree list
        graph.GetNeighbours(selected_vertex, neighbours);
        for ( int neighbour : neighbours ) {
            if ( coloring[neighbour] > 0 ) {
                continue;
            }
            list.AddNeighbourColor(neighbour, selected_color);
        }
        std::cout << "   Right after updating sat degrees: " << std::endl;
        basic_dsatur_list_test(list, 2);
        std::cout << std::endl;
    }
    
}

void test_graph(Graph& graph) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        DSaturColorStrategy color_strategy;

        unsigned short max_k;

        color_strategy.Color(graph, max_k);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();

    std::cout << "Time to color a graph with " << graph.GetNumVertices() << "vertices and " 
              << graph.GetNumEdges() << " edges: " << std::scientific << elapsed_time/std::pow(10, 9) << std::endl;

    graph.SortByColor(false);

    std::cout << "Vertices: " << TestFunctions::VecToString(graph.GetVertices()) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(graph.GetColoring()) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
    std::cout << "Is valid: " << TestFunctions::CheckColoring(graph) << std::endl;

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
    const std::string file_name = "le450_25a.col";

    std::cout << "-- TESTING DSaturList --" << std::endl;
    //test_dsatur_list("10_vertices_graph.clq");

    std::cout << "-- COLORING CSR GRAPH --" << std::endl;

    test_csr_graph("queen10_10.col");

    return 0;
}
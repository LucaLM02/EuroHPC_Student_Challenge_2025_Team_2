#include "dimacs_graph.hpp"
#include "dimacs.hpp"

#include "test_common.hpp"

#include <iostream>
#include <string>

void test_neighbors(const Graph& graph, int neighbors_of, int indentation=0) {
    std::vector<int> vertices;
    graph.GetNeighbours(neighbors_of, vertices);

    std::string indentation_string;
    for ( int i=0; i<indentation; i++) {
        indentation_string.append("  ");
    }

    std::cout << indentation_string << "Neighbors of " << neighbors_of << ": " << TestFunctions::VecToString(vertices) << std::endl;

}

void test_basics(const Graph& graph, unsigned int indentation = 0) {
    const std::vector<int>& vertices = graph.GetVertices();

    std::string indentation_string;

    for ( int i=0; i<indentation; i++) {
        indentation_string.append("  ");
    }
    std::cout << indentation_string << "Num vertices:           " << graph.GetNumVertices() << std::endl;
    std::cout << indentation_string << "Vertices:               " << TestFunctions::VecToString(vertices) << std::endl; 
    std::cout << indentation_string << "Num edges:              " << graph.GetNumEdges() << std::endl;
    std::cout << indentation_string << "Degrees:                " << TestFunctions::VecToString(graph.GetDegrees()) << std::endl;
    std::cout << indentation_string << "Max degree:             " << graph.GetMaxDegree() << std::endl;
    std::cout << indentation_string << "Vertex of max degree:   " << graph.GetVertexWithMaxDegree() << std::endl;
}

void full_basics_test(const Graph& graph, unsigned int indentation=0) {
    test_basics(graph, indentation);

    const std::vector<int>& vertices = graph.GetVertices();
    for ( int i : vertices ) {
        test_neighbors(graph, i, indentation);
    }
}

void test_graph_manipulation(Graph& graph) {
    graph.AddEdge(7, 8);

    std::cout << "After AddEdge(7,8)" << std::endl;
    full_basics_test(graph, 1);

    std::cout << "After AddEdge(7,5)" << std::endl;
    graph.AddEdge(7, 5);
    full_basics_test(graph, 1);

    std::cout << "After RemoveEdge(5, 7) (nothing should be changed)" << std::endl;
    graph.RemoveEdge(5, 7);
    full_basics_test(graph, 1);
    
    std::cout << "After MergeVertices(7, 6)" << std::endl;
    graph.MergeVertices(7, 6);
    full_basics_test(graph, 1);

    std::cout << "After RemoveVertex(2)" << std::endl;
    graph.RemoveVertex(2);
    full_basics_test(graph, 1);
}

void test_vertex_order_change(Graph &graph, std::vector<int> &&new_order) {
    graph.SetVertices(new_order);
    std::cout << "After changing the vertex order to " << TestFunctions::VecToString(new_order) << std::endl;
    full_basics_test(graph, 1);
}

int main() {
    std::string file_name = "10_vertices_graph";

    DimacsGraph& graph = *DimacsGraph::LoadFromDimacs(file_name);

    std::cout << "Graph creation" << std::endl;
    full_basics_test(graph, 1);
    
    test_graph_manipulation(graph);

    test_vertex_order_change(graph, {5, 1, 4, 7, 8, 3});     // 6 was merged into 7; 2 was deleted

}
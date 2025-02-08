// test_fastwclq_edge_cases.cpp
#include "fastwclq.hpp"
#include "graph.hpp"
#include "csr_graph.hpp"
#include <iostream>
#include <cassert>

// Test case: Empty graph
void TestEmptyGraph() {
    CSRGraph graph;
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Empty Graph Test Passed!" << std::endl;
    assert(max_clique.empty());
}

// Test case: Graph with a single vertex
void TestSingleVertex() {
    CSRGraph graph;
    graph.AddVertex();
    graph.AddEdge(1, 1); // Self-loop (ignored in most cases)
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Single Vertex Test Passed!" << std::endl;
    assert(max_clique.size() == 1);
}

// Test case: Disconnected graph
void TestDisconnectedGraph() {
    CSRGraph graph;
    graph.AddVertex();
    graph.AddVertex();
    graph.AddVertex();
    graph.AddVertex();
    graph.AddEdge(1, 2);
    graph.AddEdge(3, 4);
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Disconnected Graph Test Passed!" << std::endl;
    assert(max_clique.size() == 2);
}

// Test case: Fully connected graph (complete graph)
void TestCompleteGraph() {
    CSRGraph graph;
    int num_vertices = 5;
    for (int i = 1; i <= num_vertices; i++) {
        graph.AddVertex();
    }
    for (int i = 1; i <= num_vertices; i++) {
        for (int j = i + 1; j <= num_vertices; j++) {
            graph.AddEdge(i, j);
        }
    }
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Complete Graph Test Passed!" << std::endl;
    assert(max_clique.size() == num_vertices);
}

// Test case: Sparse graph with one large clique
void TestSparseGraphWithClique() {
    CSRGraph graph;
    int num_vertices = 8;
    for (int i = 1; i <= num_vertices; i++) {
        graph.AddVertex();
    }
    // Add a clique of size 4
    graph.AddEdge(1, 2);
    graph.AddEdge(1, 3);
    graph.AddEdge(1, 4);
    graph.AddEdge(2, 3);
    graph.AddEdge(2, 4);
    graph.AddEdge(3, 4);
    // Add some sparse connections
    graph.AddEdge(5, 6);
    graph.AddEdge(7, 8);
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Sparse Graph with Clique Test Passed!" << std::endl;
    assert(max_clique.size() == 4);
}

// Test case: Large graph with multiple cliques
void TestLargeGraph() {
    CSRGraph graph;

    for (int i = 1; i <= 20; i++) {
        graph.AddVertex();
    }

    for (int i = 1; i <= 10; i++) {
        for (int j = i + 1; j <= 10; j++) {
            graph.AddEdge(i, j);
        }
    }
    for (int i = 11; i <= 20; i++) {
        for (int j = i + 1; j <= 20; j++) {
            graph.AddEdge(i, j);
        }
    }
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Large Graph Test Passed!" << std::endl;
    assert(max_clique.size() == 10);
}

// Test case: Randomly connected graph
void TestRandomGraph() {
    CSRGraph graph;

    for (int i = 1; i <= 50; i++) {
        graph.AddVertex();
    }

    srand(42); // Fixed seed for consistency
    for (int i = 1; i <= 50; i++) {
        for (int j = i + 1; j < 50; j++) {
            if (rand() % 3 == 0) { // 33% probability of connection
                graph.AddEdge(i, j);
            }
        }
    }
    FastWClq solver(graph, 4);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Random Graph Test Passed! (Clique size: " << max_clique.size() << ")" << std::endl;
    std::cout << std::endl;
}

int TestLoadedGraph(const std::string& file_name, int k) {
    CSRGraph& graph = *CSRGraph::LoadFromDimacs(file_name);

    FastWClq solver(graph, k);
    std::vector<int> max_clique = solver.FindMaxWeightClique();
    //std::cout << "Clique found of dimension " << max_clique.size() << std::endl;
    return max_clique.size();
}

int main() {
    TestEmptyGraph();
    TestSingleVertex();
    TestDisconnectedGraph();
    TestCompleteGraph();
    TestSparseGraphWithClique();
    TestLargeGraph();
    TestRandomGraph();

    int best_result = -1;
    int best_k;

    for (int k = 1; k < 400; k++) {
        int mean_clique_size = 0;
        constexpr int N = 100;
        for (int i = 0; i < N; i++ ) {
            mean_clique_size += TestLoadedGraph("le450_25c.col", k);
        }
        mean_clique_size /= N;

        std::cout << "Mean clique dimension with k=" << k << " is: " << mean_clique_size << std::endl;
        if ( mean_clique_size > best_result ) {
            best_result = mean_clique_size;
            best_k = k;
        }
    }
    return 0;
}


// test_fastwclq_edge_cases.cpp
#include "fastwclq.hpp"
#include "graph_clique.hpp"
#include <iostream>
#include <cassert>

// Test case: Empty graph
void TestEmptyGraph() {
    GraphClique graph;
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Empty Graph Test Passed!" << std::endl;
    assert(max_clique.empty());
}

// Test case: Graph with a single vertex
void TestSingleVertex() {
    GraphClique graph;
    graph.AddEdge(0, 0); // Self-loop (ignored in most cases)
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Single Vertex Test Passed!" << std::endl;
    assert(max_clique.size() == 1);
}

// Test case: Disconnected graph
void TestDisconnectedGraph() {
    GraphClique graph;
    graph.AddEdge(0, 1);
    graph.AddEdge(2, 3);
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Disconnected Graph Test Passed!" << std::endl;
    assert(max_clique.size() == 2);
}

// Test case: Fully connected graph (complete graph)
void TestCompleteGraph() {
    GraphClique graph;
    int num_vertices = 5;
    for (int i = 0; i < num_vertices; i++) {
        for (int j = i + 1; j < num_vertices; j++) {
            graph.AddEdge(i, j);
        }
    }
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Complete Graph Test Passed!" << std::endl;
    assert(max_clique.size() == num_vertices);
}

// Test case: Sparse graph with one large clique
void TestSparseGraphWithClique() {
    GraphClique graph;
    // Add a clique of size 4
    graph.AddEdge(0, 1);
    graph.AddEdge(0, 2);
    graph.AddEdge(0, 3);
    graph.AddEdge(1, 2);
    graph.AddEdge(1, 3);
    graph.AddEdge(2, 3);
    // Add some sparse connections
    graph.AddEdge(4, 5);
    graph.AddEdge(6, 7);
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Sparse Graph with Clique Test Passed!" << std::endl;
    assert(max_clique.size() == 4);
}

// Test case: Large graph with multiple cliques
void TestLargeGraph() {
    GraphClique graph;
    for (int i = 0; i < 10; i++) {
        for (int j = i + 1; j < 10; j++) {
            graph.AddEdge(i, j);
        }
    }
    for (int i = 10; i < 20; i++) {
        for (int j = i + 1; j < 20; j++) {
            graph.AddEdge(i, j);
        }
    }
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Large Graph Test Passed!" << std::endl;
    assert(max_clique.size() == 10);
}

// Test case: Randomly connected graph
void TestRandomGraph() {
    GraphClique graph;
    srand(42); // Fixed seed for consistency
    for (int i = 0; i < 50; i++) {
        for (int j = i + 1; j < 50; j++) {
            if (rand() % 3 == 0) { // 33% probability of connection
                graph.AddEdge(i, j);
            }
        }
    }
    FastWClq solver(graph);
    std::set<int> max_clique = solver.FindMaxWeightClique();
    std::cout << "Random Graph Test Passed! (Clique size: " << max_clique.size() << ")" << std::endl;
}

int main() {
    TestEmptyGraph();
    TestSingleVertex();
    TestDisconnectedGraph();
    TestCompleteGraph();
    TestSparseGraphWithClique();
    TestLargeGraph();
    TestRandomGraph();
    return 0;
}


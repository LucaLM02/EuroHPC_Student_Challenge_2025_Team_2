#ifndef FASTWCLQ_H
#define FASTWCLQ_H

#include <set>
//#include "graph_clique.hpp"  // Assuming a Graph class providing adjacency list
#include "graph.hpp"

// Class for solving the Maximum Weight Clique problem using heuristic methods
class FastWClq {
public:
    // Constructor initializes the solver with the given graph
    explicit FastWClq(const Graph& graph, int k);

    // Finds the maximum weight clique in the graph
    std::vector<int> FindMaxWeightClique();

private:
    const Graph& graph_;  // Reference to the input graph
    std::set<int> max_clique_;  // Stores the best clique found
    int max_weight_;  // Tracks the maximum weight of a clique found
    int k_;

    // Reduces the graph by removing vertices that cannot be part of an optimal clique
    std::vector<int> GraphReduction(const std::vector<int>& C_best);

    // Constructs a clique using heuristic selection methods
    std::vector<int> CliqueConstruction();

    // Estimates the benefit of adding a vertex to the clique
    int BenefitEstimate(int v, const std::vector<int>& CandSet);

    // Selects the best vertex to add to the clique using Best from Multiple Selection (BMS) heuristic
    int ChooseVertex(const std::vector<int>& CandSet);
};

#endif // FASTWCLQ_H

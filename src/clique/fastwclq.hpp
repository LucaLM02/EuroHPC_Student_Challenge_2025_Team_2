#ifndef FASTWCLQ_H
#define FASTWCLQ_H

#include <set>
#include "graph_clique.hpp"  // Assuming a Graph class providing adjacency list

// Class for solving the Maximum Weight Clique problem using heuristic methods
class FastWClq {
public:
    // Constructor initializes the solver with the given graph
    explicit FastWClq(const GraphClique& graph);

    // Finds the maximum weight clique in the graph
    std::set<int> FindMaxWeightClique();

private:
    const GraphClique& graph_;  // Reference to the input graph
    std::set<int> max_clique_;  // Stores the best clique found
    int max_weight_;  // Tracks the maximum weight of a clique found

    // Reduces the graph by removing vertices that cannot be part of an optimal clique
    std::set<int> GraphReduction(const std::set<int>& C_best);

    // Constructs a clique using heuristic selection methods
    std::set<int> CliqueConstruction();

    // Estimates the benefit of adding a vertex to the clique
    int BenefitEstimate(int v, const std::set<int>& CandSet);

    // Selects the best vertex to add to the clique using Best from Multiple Selection (BMS) heuristic
    int ChooseVertex(const std::set<int>& CandSet, int k);
};

#endif // FASTWCLQ_H

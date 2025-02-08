#include "fastwclq.hpp"
#include <algorithm>
#include <iterator>
#include <random>
#include <vector>
#include <iostream>
#include <set>
#include <unordered_set>
#include <unordered_map>

// Constructor initializes the graph reference and sets the max weight to zero
FastWClq::FastWClq(const GraphClique& graph) : graph_(graph), max_weight_(0) {}

// Main function to find the maximum weight clique
std::set<int> FastWClq::FindMaxWeightClique() {
    std::set<int> best_clique;
    int iteration = 0;

    while (true) {
        iteration++;

        // Construct a clique using heuristic methods
        std::set<int> clique = CliqueConstruction();
        std::cout << "Clique Size Found: " << clique.size() << std::endl;

        // Update the best found clique if the new one is larger
        if (clique.size() > best_clique.size()) {
            best_clique = clique;
        }

        // Reduce the graph based on the best clique found so far
        std::set<int> reduced_graph = GraphReduction(best_clique);
        std::cout << "Reduced Graph Size: " << reduced_graph.size() << std::endl;

        // Check if the reduced graph is empty
        if (reduced_graph.empty() || reduced_graph.size() == 1) {
            std::cout << "Terminating: Final Max Clique Size = " << best_clique.size() << std::endl;
            return best_clique;
        }



        // Prevent infinite loops (temporary safeguard)
        if (iteration > 100) {
            std::cerr << "ERROR: Algorithm stuck in an infinite loop. Exiting." << std::endl;
            exit(1);
        }
    }
}


// Graph reduction step: removes vertices that cannot be part of a maximum clique
std::set<int> FastWClq::GraphReduction(const std::set<int>& C_best) {
    std::set<int> reduced_graph;
    int current_best_size = C_best.size();

    for (int v : graph_.GetVertices()) {
        int upper_bound = graph_.GetVertexWeight(v) + graph_.GetNeighborsWeightSum(v);

        // Ensure only valid vertices remain in the reduced graph
        if (upper_bound > current_best_size) {
            reduced_graph.insert(v);
        }

        // Ensure algorithm terminates when only one node remains
        if (reduced_graph.size() <= 1) {
            std::cout << "Graph reduced to a single node. Terminating." << std::endl;
            return {};
        }

    }

    // Prevent the same graph from being passed without reduction
    if (reduced_graph == graph_.GetVertices()) {
        std::cerr << "Warning: No vertices were removed in graph reduction." << std::endl;
        return {};
    }

    return reduced_graph;
}


// Construct a clique iteratively using a greedy selection method
std::set<int> FastWClq::CliqueConstruction() {
    std::set<int> C;
    std::set<int> CandSet = graph_.GetVertices();

    while (!CandSet.empty()) {
        int v = ChooseVertex(CandSet, 4);
        C.insert(v);

        // Handle case where there's only one vertex in the graph
        if (CandSet.size() == 1) {
            break;
        }

        std::set<int> new_CandSet;
        for (int u : CandSet) {
            if (graph_.IsNeighbor(v, u)) {
                new_CandSet.insert(u);
            }
        }
        CandSet = new_CandSet;
    }

    return C;
}

// Estimate the benefit of adding vertex v to the clique
int FastWClq::BenefitEstimate(int v, const std::set<int>& CandSet) {
    return graph_.GetVertexWeight(v) + graph_.GetNeighborsWeightSum(v) / 2;
}

// Selects the best vertex to add to the clique using a heuristic method
int FastWClq::ChooseVertex(const std::set<int>& CandSet, int k) {
    // If the candidate set is small, choose the best directly
    if (CandSet.size() <= k) {
        return *std::max_element(CandSet.begin(), CandSet.end(), [&](int a, int b) {
            return BenefitEstimate(a, CandSet) < BenefitEstimate(b, CandSet);
        });
    }

    // Use random sampling - BMS heuristic - to pick a subset of k candidates
    std::vector<int> samples;
    std::sample(CandSet.begin(), CandSet.end(), std::back_inserter(samples), k, std::mt19937{std::random_device{}()});

    // Return the best among the sampled candidates based on benefit estimate
    return *std::max_element(samples.begin(), samples.end(), [&](int a, int b) {
        return BenefitEstimate(a, CandSet) < BenefitEstimate(b, CandSet);
    });
}

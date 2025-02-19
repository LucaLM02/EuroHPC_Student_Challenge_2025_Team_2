#include "fastwclq.hpp"
#include <algorithm>
#include <iterator>
#include <random>
#include <vector>
#include <iostream>
#include <set>
#include <unordered_set>
#include <unordered_map>


int FastCliqueStrategy::FindClique(const Graph &graph) const
{

    _solver = std::make_unique<FastWClq>(graph, _k);
    
    return _solver->FindMaxWeightClique().size();
}

std::vector<int> FastCliqueStrategy::GetClique() const
{
    return _solver->GetMaxClique();
}

// Constructor initializes the graph reference and sets the max weight to zero
FastWClq::FastWClq(const Graph& graph, int k) : graph_(graph), max_weight_(0), k_{k} {}

// Main function to find the maximum weight clique
std::vector<int> FastWClq::FindMaxWeightClique() {
    std::vector<int> best_clique;
    int iteration = 0;
    
    while (true) {
        iteration++;
        // Construct a clique using heuristic methods
        std::vector<int> clique = CliqueConstruction();
        // Update the best found clique if the new one is larger
        if (clique.size() > best_clique.size()) {
            best_clique = clique;
        }

        // Reduce the graph based on the best clique found so far
        std::vector<int> reduced_graph = GraphReduction(best_clique);
        // Check if the reduced graph is empty
        if (reduced_graph.empty() || reduced_graph.size() == 1) {
            max_clique_ = best_clique;
            return best_clique;
        }

        // Prevent infinite loops (temporary safeguard)
        if (iteration > 100) {
            //std::cerr << "ERROR: Algorithm stuck in an infinite loop. Exiting." << std::endl;
            max_clique_ = best_clique;
            return best_clique;
        }
    }
}


// Graph reduction step: removes vertices that cannot be part of a maximum clique
std::vector<int> FastWClq::GraphReduction(const std::vector<int>& C_best) {
    //std::set<int> reduced_graph;
    std::vector<int> reduced_graph;
    int current_best_size = C_best.size();


    for (int v : graph_.GetVertices()) {
        int upper_bound = 1 + graph_.GetExDegree(v);

        // Ensure only valid vertices remain in the reduced graph
        if (upper_bound > current_best_size) {
            //reduced_graph.insert(v);
            reduced_graph.push_back(v);
        }

        // Ensure algorithm terminates when only one node remains
        if (reduced_graph.size() <= 1) {
            return {};
        }

    }

    // Prevent the same graph from being passed without reduction
    //if (reduced_graph == graph_.GetVertices()) {
    if (reduced_graph.size() == graph_.GetVertices().size()) {
        //std::cerr << "Warning: No vertices were removed in graph reduction." << std::endl;
        return {};
    }

    return reduced_graph;
}


// Construct a clique iteratively using a greedy selection method
//std::set<int> FastWClq::CliqueConstruction() {
std::vector<int> FastWClq::CliqueConstruction() {
    //std::set<int> C;
    //std::set<int> CandSet = graph_.GetVertices();
    std::vector<int> C;     // size?
    std::vector<int> CandSet = graph_.GetVertices();

    while (!CandSet.empty()) {
        int v = ChooseVertex(CandSet);
        C.push_back(v);

        // Handle case where there's only one vertex in the graph
        if (CandSet.size() == 1) {
            break;
        }

        std::vector<int> new_CandSet;
        new_CandSet.reserve(CandSet.size());
        for (int u : CandSet) {
            //if (graph_.IsNeighbor(v, u)) {
                if (v < 0 || u < 0) {
                    std::cerr << "error invalid vertex" << v << " u=" << u << std::endl;
                    continue;
                }
            if (graph_.HasEdge(v, u)) {
                new_CandSet.push_back(u);
            }
        }
        CandSet = new_CandSet;
    }

    return C;
}

// Estimate the benefit of adding vertex v to the clique
int FastWClq::BenefitEstimate(int v, const std::vector<int>& CandSet) {
    //return graph_.GetVertexWeight(v) + graph_.GetNeighborsWeightSum(v) / 2;
    return graph_.GetDegree(v) + graph_.GetExDegree(v) / 2;
}

// Selects the best vertex to add to the clique using a heuristic method
//int FastWClq::ChooseVertex(const std::set<int>& CandSet, int k) {
int FastWClq::ChooseVertex(const std::vector<int>& CandSet) {
    // If the candidate set is small, choose the best directly
    if (CandSet.size() <= k_) {
        int result = *std::max_element(CandSet.begin(), CandSet.end(), [&](int a, int b) {
            return BenefitEstimate(a, CandSet) < BenefitEstimate(b, CandSet);
        });
        return result;
    }

    // Use random sampling - BMS heuristic - to pick a subset of k candidates
    std::vector<int> samples;
    std::sample(CandSet.begin(), CandSet.end(), std::back_inserter(samples), k_, std::mt19937{std::random_device{}()});

    // Return the best among the sampled candidates based on benefit estimate
    
    int result = *std::max_element(samples.begin(), samples.end(), [&](int a, int b) {
        return BenefitEstimate(a, CandSet) < BenefitEstimate(b, CandSet);
    });
    return result;
}

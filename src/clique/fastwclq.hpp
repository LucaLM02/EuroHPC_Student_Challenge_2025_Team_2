#ifndef FASTWCLQ_HPP
#define FASTWCLQ_HPP

#include <set>
#include <memory>

#include "graph.hpp"
#include "clique_strategy.hpp"

class FastWClq;

class FastCliqueStrategy : public CliqueStrategy {
    public:
        FastCliqueStrategy(int k=5) : _k{k} {};
        /**
         * @brief finds the size of a feasible clique in the graph.
         * 
         * @param graph the graph of which the clique has to be found
         * @returns the size of the feasible clique
         */
        virtual int FindClique(const Graph &graph) const override;
        virtual std::vector<int> GetClique() const override;
    private:
        mutable std::unique_ptr<FastWClq> _solver;
        const int _k;
};

// Class for solving the Maximum Weight Clique problem using heuristic methods
class FastWClq {
public:
    // Constructor initializes the solver with the given graph
    explicit FastWClq(const Graph& graph, int k);

    // Finds the maximum weight clique in the graph
    std::vector<int> FindMaxWeightClique();

    std::vector<int> GetMaxClique() { return max_clique_; }

private:
    const Graph& graph_;  // Reference to the input graph
    std::vector<int> max_clique_;  // Stores the best clique found
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

#endif // FASTWCLQ_HPP

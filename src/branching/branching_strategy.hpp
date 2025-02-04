#ifndef BRANCHING_STRATEGY_HPP
#define BRANCHING_STRATEGY_HPP

#include "common.hpp"
#include "graph.hpp"

#include <random>
#include <memory>

class BranchingStrategy {
    public:
        BranchingStrategy() 
        {}

        virtual std::pair<unsigned int, unsigned int> ChooseVertices(const Graph& graph) = 0;
};

class RandomBranchingStrategy : public BranchingStrategy {
    public:
        RandomBranchingStrategy(int num_vertices);

        virtual std::pair<unsigned int, unsigned int> ChooseVertices(const Graph& graph) override;

    protected:
        std::pair<unsigned int, unsigned int> _vertex_pair;
        std::unique_ptr<std::mt19937> _random_generator;

    private:
};

class DegreeBranchingStrategy : public BranchingStrategy {
    public:
        DegreeBranchingStrategy();

        virtual std::pair<unsigned int, unsigned int> ChooseVertices(const Graph& graph) override;

};

/**
 * @brief prefers choosing the elements of the clique, by increasing the clique size
 * 
 */
class CliqueBranchingStrategy : public BranchingStrategy {
    public:
        CliqueBranchingStrategy(/*const CliqueStrategy& clique*/);

        virtual std::pair<unsigned int, unsigned int> ChooseVertices(const Graph& graph) override;
};

#endif // BRANCHING_STRATEGY_HPP

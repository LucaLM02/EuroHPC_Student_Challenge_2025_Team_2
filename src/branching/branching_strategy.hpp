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

        virtual std::pair<unsigned int, unsigned int> ChooseVertices(const Graph& graph);
    protected:
        std::pair<unsigned int, unsigned int> _vertex_pair;

        virtual void _PermuteVertices(const Graph& graph) = 0;
};

class RandomBranchingStrategy : public BranchingStrategy {
    public:
        RandomBranchingStrategy(int num_vertices) 
        : BranchingStrategy(), _uniform_distribution(0, num_vertices-1)
        {
            std::random_device dev;
            _random_generator = std::make_unique<std::mt19937>(dev());
        }
    
    protected:
        virtual void _PermuteVertices(const Graph& graph) override;

        std::unique_ptr<std::mt19937> _random_generator;

    private:
        std::uniform_int_distribution<std::mt19937::result_type> _uniform_distribution;
};

#endif // BRANCHING_STRATEGY_HPP

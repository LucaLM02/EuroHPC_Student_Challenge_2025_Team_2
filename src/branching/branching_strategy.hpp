#ifndef BRANCHING_STRATEGY_HPP
#define BRANCHING_STRATEGY_HPP

#include "common.hpp"
#include "graph.hpp"

#include <random>
#include <memory>

/**
 * @brief abstact strategy with which, at each step of the branch and bound, a new 
 *        pair of vertices to merge (contract) or to set as neighbours is chosen.
 *        It is a wrapper of method ChooseVertices.
 *        It is an abstract class, concrete strategies must inherit from this
 * 
 */
class BranchingStrategy {
    public:
        
        BranchingStrategy() 
        {}

        virtual std::pair<int, int> 
        ChooseVertices(const Graph& graph) = 0;
};

class RandomBranchingStrategy : public BranchingStrategy {
    public:
        RandomBranchingStrategy(int num_vertices);

        virtual std::pair<int, int> 
        ChooseVertices(const Graph& graph) override;

    protected:
        std::pair<unsigned int, unsigned int> _vertex_pair;
        std::unique_ptr<std::mt19937> _random_generator;

    private:
};

class DegreeBranchingStrategy : public BranchingStrategy {
    public:
        DegreeBranchingStrategy();

        virtual std::pair<int, int> 
        ChooseVertices(const Graph& graph) override;

};

/**
 * @brief the vertices are selected from certain maximal sets which are determined at the beginning.
 * @details the selection of vertices can be customized in the following ways:
 *      - number of times elements of the same sets are chosen (yelds Type::Equal pair)
 *      - number of times elements of the different sets are chosen (yelds Type::Diff pair)
 *      - how frequently updating the sets (and if not uupdate sets)
 *      - how many sets to find
 * 
 */
class IndependentSetBranchingStrategy : public BranchingStrategy {
    public:
        IndependentSetBranchingStrategy();

        /**
         * @brief sets the frequency for recalculating the independent sets
         * @details since the independent sets get smaller, after certain time it might be appropriate to
         *          find new independent sets
         * 
         * @param frequency 
         */
        void SetUpdateFrequency(unsigned int frequency);

        /**
         * @brief sets how many consequent times two elements of different sets are selected
         * 
         * @param length sequence length 
         */
        void SetLengthDiffSequence(unsigned int length);

        /**
         * @brief sets how many consequent times two elements of the same set are selected
         * 
         * @param length sequence length
         */
        void SetLengthEqualSequence(unsigned int length);

        /**
         * @brief finds the maximal independent sets of the graph
         * @param graph graph where to search the sets
         * @param number number of sets to search
         */
        void FindIndependentSets(const Graph& graph, unsigned int number);
        virtual std::pair<int, int> 
        ChooseVertices(const Graph& graph) override;
    
    protected:

        unsigned int _iteration_counter;
        unsigned int _update_frequency;

        bool _is_equal_sequence;
        unsigned int _current_length;
        unsigned int _length_equal;
        unsigned int _length_diff;

        int _chosen_indep_set;

    private: 
        std::vector<std::vector<int>> _independent_sets;
};

/**
 * @brief prefers choosing the elements of the clique, by increasing the clique size
 * 
 */
class CliqueBranchingStrategy : public BranchingStrategy {
    public:
        CliqueBranchingStrategy(/*const CliqueStrategy& clique*/);

        virtual std::pair<int, int> 
        ChooseVertices(const Graph& graph) override;
};

/**
 * @brief strategy with which, at each step of the branch and bound, a new 
 *        pair of vertices to merge (contract) or to set as neighbours is chosen.
 *        In particular, vertices x and y are choosen such that they are the pair
 *        of non-adjacent vertices with the highest number of common neighbours
 * 
 */
class NeighboursBranchingStrategy : public BranchingStrategy {
    public:
        NeighboursBranchingStrategy() = default;

        virtual std::pair<int, int> 
        ChooseVertices(const Graph& graph) override;
};

#endif // BRANCHING_STRATEGY_HPP

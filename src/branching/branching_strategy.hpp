#ifndef BRANCHING_STRATEGY_HPP
#define BRANCHING_STRATEGY_HPP

#include "common.hpp"
#include "graph.hpp"

#include <random>
#include <memory>

class BranchingStrategy {
    public:
        enum PairType {
            Equal,
            Diff,
            DontCare
        };
        
        BranchingStrategy() 
        {}

        virtual std::pair<unsigned int, unsigned int> 
        ChooseVertices(const Graph& graph, PairType& type) = 0;
};

class RandomBranchingStrategy : public BranchingStrategy {
    public:
        RandomBranchingStrategy(int num_vertices);

        virtual std::pair<unsigned int, unsigned int> 
        ChooseVertices(const Graph& graph, PairType& type) override;

    protected:
        std::pair<unsigned int, unsigned int> _vertex_pair;
        std::unique_ptr<std::mt19937> _random_generator;

    private:
};

class DegreeBranchingStrategy : public BranchingStrategy {
    public:
        DegreeBranchingStrategy();

        virtual std::pair<unsigned int, unsigned int> 
        ChooseVertices(const Graph& graph, PairType& type) override;

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
        virtual std::pair<unsigned int, unsigned int> 
        ChooseVertices(const Graph& graph, PairType& type) override;
    
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

        virtual std::pair<unsigned int, unsigned int> 
        ChooseVertices(const Graph& graph, PairType& type) override;
};

#endif // BRANCHING_STRATEGY_HPP

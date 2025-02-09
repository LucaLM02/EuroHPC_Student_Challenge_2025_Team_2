#ifndef CLIQUESTRATEGY_HPP
#define CLIQUESTRATEGY_HPP

#include "graph.hpp"

/**
 *  @brief wrapper class for the FindClique methods
 *  @see the description of methods
 */
class CliqueStrategy {
    public:
        /**
         * @brief finds the size of a feasible clique in the graph.
         * 
         * @param graph the graph of which the clique has to be found
         * @returns the size of the feasible clique
         */
        virtual int FindClique(const Graph &graph) const = 0;
        /**
         * @brief returns the last computed maximal clique
         * 
         * @returns the last computed maximal clique
         */
        virtual std::vector<int> GetClique() const = 0;
};

/**
 * @brief temporary class used only for testing CliqueStrategy
 * 
 */
class StubCliqueStrategy : public CliqueStrategy {
    public:
        virtual int FindClique(const Graph &graph) const override;
        virtual std::vector<int> GetClique() const override { return {}; }
    
};

#endif // CLIQUESTRATEGY_HPP

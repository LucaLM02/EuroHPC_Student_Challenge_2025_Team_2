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
     * @brief finds a clique of the graph.
     * 
     * @param graph the graph of which the clique has to be found
     * @param clique the vertices belonging to the clique
     */
    virtual void FindClique(const Graph &graph, 
                            std::vector<int>& clique) const = 0;

    /**
     * @brief finds the size of a feasible clique in the graph.
     * 
     * @param graph the graph of which the clique has to be found
     * @returns the size of the feasible clique
     */
    virtual int FindClique(const Graph &graph) const = 0;
};

/**
 * @brief temporary class used only for testing CliqueStrategy
 * 
 */
class StubCliqueStrategy : public CliqueStrategy {
    public:
    void FindClique(const Graph& graph,
                    std::vector<int>& clique) const;
    
    virtual int FindClique(const Graph &graph) const;
};

#endif // CLIQUESTRATEGY_HPP

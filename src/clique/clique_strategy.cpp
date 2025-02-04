#include "clique_strategy.hpp"

void StubCliqueStrategy::FindClique(const Graph &graph, std::vector<int> &clique) const
{

}

int StubCliqueStrategy::FindClique(const Graph &graph) const
{
    if ( graph.GetNumEdges() >= 1 ) {
        return 2;
    }
    return 1;
}

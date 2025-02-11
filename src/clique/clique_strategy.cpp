#include "clique_strategy.hpp"

int StubCliqueStrategy::FindClique(const Graph &graph) const
{
    if ( graph.GetNumEdges() >= 1 ) {
        return 2;
    }
    return 1;
}

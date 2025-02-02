#include "branching_strategy.hpp"

std::pair<unsigned int, unsigned int> BranchingStrategy::ChooseVertices(const Graph& graph) {
    _PermuteVertices(graph);

    return _vertex_pair;
}

void RandomBranchingStrategy::_PermuteVertices(const Graph& graph) {
    do 
    {
        _vertex_pair.first   = _uniform_distribution(*_random_generator);
        _vertex_pair.second  = _uniform_distribution(*_random_generator);

    } while (graph.HasEdge(_vertex_pair.first, _vertex_pair.second));
}

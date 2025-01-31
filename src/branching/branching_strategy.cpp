#include "branching_strategy.hpp"

std::pair<unsigned int, unsigned int> BranchingStrategy::ChooseVertices(Constraints& config) {
    _PermuteVertices(config);

    return _vertex_pair;
}

void RandomBranchingStrategy::_PermuteVertices(Constraints& config) {
    do 
    {
        _vertex_pair.first   = _uniform_distribution(*_random_generator);
        _vertex_pair.second  = _uniform_distribution(*_random_generator);

    } while (config.IsEqualConstraint(_vertex_pair) || config.IsDiffConstraint(_vertex_pair));
}

#include "minimum_chromatic_finder.hpp"


// lower_bound and upper_bound should be global variables
void MinimumChromaticFinder::Solve() {

    std::vector<unsigned short> coloring;
    unsigned short k_max;
    //_color_strategy.Color(_vertices, _edges, _current_config, coloring, k_max);
    _color_strategy.Color(_vertices, _edges, _current_config, coloring, k_max, _vertices.size());

    // tries to reduce the color of the worst vertices
    // might not be the best strategy
    _recolor_strategy.Recolor(_vertices, _edges, _current_config, coloring);

    if ( k_max < _upper_bound ) {
        _upper_bound = k_max;
        _best_config = _current_config;
    }

    // what about updating the lower bound??

    std::pair<unsigned int, unsigned int> vertex_pair   = _branching_strategy.ChooseVertices(_current_config);

    // 1st branch: same color
    _current_config.AddEqualConstraint(vertex_pair);
    Solve();

    // 2nd branch: different color
    _current_config.AddDiffConstraint(vertex_pair);
    Solve();

}
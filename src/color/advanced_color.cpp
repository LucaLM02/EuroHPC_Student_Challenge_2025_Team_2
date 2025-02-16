#include "advanced_color.hpp"

void InterleavedColorStrategy::Color(Graph &graph, unsigned short &k_max) const
{
    _curr_length++;
    if ( _is_first ) {
        _first_color_strategy.Color(graph, k_max);
        if ( _curr_length == _length_first_sequence ) {
            _is_first = false;
            _curr_length = 0;
        }
    } else {
        _second_color_strategy.Color(graph, k_max);
        if ( _curr_length == _length_second_sequence ) {
            _is_first = true;
            _curr_length = 0;
        }
    }
}

void InactiveColorStrategy::Color(Graph &graph, unsigned short &k_max) const
{
    std::vector<unsigned short> coloring = graph.GetColoring();

    k_max = *std::max_element(coloring.begin(), coloring.end());
}

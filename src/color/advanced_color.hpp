#ifndef ADVANCED_COLOR_HPP
#define ADVANCED_COLOR_HPP

#include "color.hpp"
#include "recolor.hpp"

class InactiveColorStrategy : public ColorStrategy {
        void Color(Graph &graph,
                    unsigned short& k_max) const;
};

class InterleavedColorStrategy : public ColorStrategy {
    public:
        InterleavedColorStrategy(ColorStrategy& first_color_strategy,
                                 ColorStrategy& second_color_strategy,
                                 unsigned int length_first_sequence, 
                                 unsigned int length_second_sequence)
            : _first_color_strategy{first_color_strategy},
              _second_color_strategy{second_color_strategy},
              _length_first_sequence{length_first_sequence},
              _length_second_sequence{length_second_sequence}
        {}

        void Color(Graph &graph,
                    unsigned short& k_max) const override {
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

    private:
        mutable bool _is_first = true;
        ColorStrategy& _first_color_strategy;
        ColorStrategy& _second_color_strategy;
        const unsigned int _length_first_sequence;
        const unsigned int _length_second_sequence;
        mutable unsigned int _curr_length;
};

#endif // ADVANCED_COLOR_HPP

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
                                 unsigned int length_color_sequence, 
                                 unsigned int length_inactive_sequence)
            : _first_color_strategy{first_color_strategy},
              _second_color_strategy{second_color_strategy},
              _length_first_sequence{length_color_sequence},
              _length_second_sequence{length_inactive_sequence}
        {}

        void Color(Graph &graph,
                    unsigned short& k_max) const;

    private:
        mutable bool _is_first = true;
        ColorStrategy& _first_color_strategy;
        ColorStrategy& _second_color_strategy;
        const unsigned int _length_first_sequence;
        const unsigned int _length_second_sequence;
        mutable unsigned int _curr_length;
};

#endif // ADVANCED_COLOR_HPP

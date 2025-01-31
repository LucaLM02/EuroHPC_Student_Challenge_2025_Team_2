#ifndef MINIMUM_CHROMATIC_FINDER_HPP
#define MINIMUM_CHROMATIC_FINDER_HPP

#include "common.hpp"
#include "color.hpp"
#include "recolor.hpp"
#include "branching_strategy.hpp"

class MinimumChromaticFinder {
    public:
        MinimumChromaticFinder(const VertexSet& vertices, const Edges& edges, 
                              BranchingStrategy& branching_strategy,
                              const ColorStrategy& color_strategy,
                              const RecolorStrategy& recolor_strategy)
        : _vertices{vertices}, _edges{edges}, 
          _branching_strategy{branching_strategy},
          _color_strategy{color_strategy},
          _recolor_strategy{recolor_strategy},
          _best_config(vertices.size())
        {}

        void Solve();
    private:
        VertexSet              _vertices;
        const Edges&           _edges;
        BranchingStrategy&     _branching_strategy;
        const ColorStrategy&   _color_strategy;
        const RecolorStrategy& _recolor_strategy;

        unsigned int _lower_bound;
        unsigned int _upper_bound;

        Constraints _best_config;
        Constraints _current_config;
};


#endif // MINIMUM_CHROMATIC_FINDER_HPP

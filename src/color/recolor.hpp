#ifndef RECOLOR_HPP
#define RECOLOR_HPP

#include "common.hpp"

/*
    @brief abstract functional class that wraps Recolor method; (partially) recolors a graph 
           to reduce the maximum color number
    See the method description
*/
class RecolorStrategy {
    public:
        /*
            @brief tries to reduce the maximum color number of a certain coloring

            @param vertices     set of vertices which compose the graph
            @param edges        set of edges which compose the graph
            @param coloring     full coloring of the graph

            @returns            0 if recoloring failed, otherwise how much the highest color was reduced
        */
        virtual unsigned int Recolor(VertexSet& vertices, const Edges& edges, 
                                    std::vector<unsigned short>& coloring) const = 0;
};

/*
    @brief functional class for the Recolor method; (partially) recolors a graph to reduce the maximum color number
    See the method description
*/
class GreedySwapRecolorStrategy : public RecolorStrategy {
    public:
        /*
            @brief tries to reduce the maximum color number of a certain coloring by greedly 
                swapping colors

            @details For each vertex with the highest degree, reduces the color by swapping 
                     it with a neighbour, which will be given another color. 
                     When all the vertices with the highest degree are reduced to a lower degree,
                     a local variable, `reduction_counter`, is incremented by one and the algorithm
                     proceeds with the new highest degree.
                     If any of this passages is not feasible, then algorithm stops and the current 
                     value of `reduction_counter` is returned

            @param vertices     set of vertices which compose the graph
            @param edges        set of edges which compose the graph
            @param coloring     full coloring of the graph

            @returns            0 if recoloring failed, otherwise how much the highest color was reduced
        */
        unsigned int Recolor(VertexSet& vertices, const Edges& edges, 
                            std::vector<unsigned short>& coloring) const override;
};

#endif // RECOLOR_HPP

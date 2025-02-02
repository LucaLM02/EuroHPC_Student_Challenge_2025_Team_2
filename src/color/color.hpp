#ifndef COLOR_HPP
#define COLOR_HPP

#include "common.hpp"

/*
    @brief functional class that wraps Color method. Colors a graph in such a way that no 
           vertices have the same color
    See the method description
*/
class ColorStrategy {
    public:
        /*
            @details Colors a graph in such a way that no vertices have the same color
                     Colors are contiguosly used from k=1 to `k_max`, where `k_max` is 
                     decided dynamically by this method

            @param vertices               vertices of the graph. Can be modified by the method
            @param edges                  edges of the graph
            @param coloring               the color for each vertex, using the same order of vertices
            @param k_mx                   highest color used
            @param expected_maximum_color expected maximum color
        */
        virtual void Color(VertexSet& vertices, const Edges& edges, 
                           std::vector<unsigned short>& coloring, 
                           unsigned short& k_max,
                           const unsigned int& expected_maximum_color) const = 0;
};

/*
    @brief functional class that wraps Color method. Colors a graph in such a way that no 
           vertices have the same color. Color is choosen in a greedy fashion, i.e. for 
           each node the lowest color is chosen
    See the method description
*/
class GreedyColorStrategy : public ColorStrategy {
    public:
        /*
            @details Colors a graph in such a way that no vertices have the same color
                     Colors are contiguosly used from k=1 to `k_max`, where `k_max` is 
                     decided dynamically by this method
                     For each vertex the lowest possible color is chosen
            

            @param vertices               vertices of the graph. Can be modified by the method
            @param edges                  edges of the graph
            @param constraints            constraints on the color that a vertex of the given graph can assume
            @param coloring               the color for each vertex, using the same order of vertices
            @param k_mx                   highest color used
            @param expected_maximum_color expected maximum color
        */
        void Color(VertexSet& vertices, const Edges& edges, 
                   std::vector<unsigned short>& coloring, 
                   unsigned short& max_k,
                   const unsigned int& expected_maximum_color) const override;
};


#endif // COLOR_HPP

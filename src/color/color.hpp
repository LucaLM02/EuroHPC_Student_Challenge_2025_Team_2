#ifndef COLOR_HPP
#define COLOR_HPP

#include "common.hpp"
#include "graph.hpp"

/**
 *  @brief functional class that wraps Color method. Colors a graph in such a way that no 
 *         vertices have the same color
 *  @see the description of method Color
 */
class ColorStrategy {
    public:
        /**
          * @details Colors a graph in such a way that no vertices have the same color
          *          Colors are contiguosly used from k=1 to `k_max`, where `k_max` is 
          *          decided dynamically by this method
          * 
          *  @param graph        the graph to color
          *  @param coloring     the color for each vertex, using the same order of vertices
          *  @param k_mx         highest color used
          */
        virtual void Color(Graph &graph,
                           std::vector<unsigned short>& coloring, 
                           unsigned short& k_max) const = 0;
};

/**
 *  @brief functional class that wraps Color method. Colors a graph in such a way that no 
 *         vertices have the same color. Color is choosen in a greedy fashion, i.e. for 
 *         each node the lowest color is chosen
 *  @see the description of method Color
 */
class GreedyColorStrategy : public ColorStrategy {
    public:
        /**
         *  @details Colors a graph in such a way that no vertices have the same color
         *           Colors are contiguosly used from k=1 to `k_max`, where `k_max` is 
         *           decided dynamically by this method
         *           For each vertex the lowest possible color is chosen
         *  
         *  @param graph        the graph to color
         *  @param coloring     the color for each vertex, using the same order of vertices
         *  @param k_mx         highest color used
         */
        void Color(Graph& graph,
                   std::vector<unsigned short>& coloring, 
                   unsigned short& max_k) const override;
};


#endif // COLOR_HPP

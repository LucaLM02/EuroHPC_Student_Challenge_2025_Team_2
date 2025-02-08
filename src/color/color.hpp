#ifndef COLOR_HPP
#define COLOR_HPP

#include "common.hpp"
#include "graph.hpp"

/**
 *  @brief functional class that wraps Color method. Colors a graph in such a way that 
 *         no adjacent vertices have the same color
 *  @see the description of method Color
 */
class ColorStrategy {
    public:
        /**
         * @details Colors a graph in such a way that no adjacent vertices have the 
         *          same color
         *          Colors are contiguosly used from k=1 to `k_max`, where `k_max` is 
         *          decided dynamically by this method
         *          Coloring is set to the graph and accessible with `Graph::GetColoring()`
         * 
         *  @param graph        the graph to color
         *  @param k_mx         highest color used
         */
        virtual void Color(Graph &graph,
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
         *           Coloring is set to the graph and accessible with `Graph::GetColoring()`
         *  
         *  @param graph        the graph to color
         *  @param k_mx         highest color used
         */
        void Color(Graph& graph,
                   unsigned short& max_k) const override;
};

/**
 * @brief finds the lowest color available for a vertex given the current (partial) coloring
 * 
 * @param graph graph which is being colored
 * @param vertex vertex which has to be colored
 * @param coloring current partial coloring
 * @param current_max_k maximum color to choose from
 * @return unsigned int color choosen
 */
unsigned short GreedyFindColor(const Graph& graph,
                             const unsigned int vertex,
                             std::vector<unsigned short>& coloring, 
                             unsigned int current_max_k);


#endif // COLOR_HPP

#ifndef DSATUR_COLOR_HPP
#define DSATUR_COLOR_HPP

#include <list>

#include "color.hpp"

/**
 *  @brief functional class that wraps Color method. Colors a graph in such a way that no 
 *         vertices have the same color. Color is choosen in a greedy fashion, but the 
 *         next vertex is the one with highest saturation degree
 *  @see the description of method Color
 */
class DSaturColorStrategy : public ColorStrategy {
    public:
        /**
         *  @details Colors a graph in such a way that no vertices have the same color
         *           Colors are contiguosly used from k=1 to `k_max`, where `k_max` is 
         *           decided dynamically by this method <br>
         *           For each vertex the lowest possible color is chosen <br>
         *           Vertex are dynamically choosen; the one with highest saturation
         *           degree is choosen. The saturation degree is the number of 
         *           different colors that the neighbors have
         *           Coloring is set to the graph and accessible with `Graph::GetColoring()`
         *  
         *  @param graph        the graph to color
         *  @param k_mx         highest color used
         */
        void Color(Graph& graph,
                   unsigned short& max_k) const override;
};

/**
 * @brief an item of the linked list used in the DSaturList
 */
typedef
struct _DSaturItem {
    _DSaturItem(int sat_degree_, int degree_, int vertex_, _DSaturItem* prev_)
    : sat_degree{sat_degree_}, degree{degree_}, vertex{vertex_}, prev{prev_}
    {}

    int sat_degree;
    int degree;
    int vertex;
    _DSaturItem* next;
    _DSaturItem* prev;
} DSaturItem;

/**
 * @brief data structure to optimize the time complexity of DSaturColorStrategy
 * @warning memory complexity is very bad since for each vertex it is stored a pointer 
 *          to a DSaturItem and this contains 2 other pointers, both degrees, the vertex
 *          and all its neighbour colors using a set, which is memory consuming
 * @note opposite to what written in the warning, note that even though the memory used
 *       is bad with respect to the number of vertices, in a graph the number of vertices
 *       is, usually, one order smaller than the number of edges. 
 *       This means that it has a limited impact with respect to a Graph class
 */
class DSaturList {
    public:
        DSaturList(const Graph& graph);
        virtual ~DSaturList();
        /**
         * @brief adds a color to the neighbour color set of that particular vertex.
         *        If it's saturation degree increases, then the underlying data structure
         *        is updated accordingly
         */
        void AddNeighbourColor(int vertex, unsigned short color);
        /**
         * @brief rerturns the saturation degree corresponding to the GetLowestVertex()
         * 
         * @returns the saturation degree corresponding to the GetLowestVertex()
         */
        int GetLowestSatDegree() const;
        /**
         * @brief returns the vertex with the lowest degree and lowest saturation degree
         * 
         * @return the vertex
         */
        int GetLowestVertex() const;
        /**
         * @brief returns and remove the vertex with the lowest degree and lowest 
         *        saturation degree
         * 
         * @return the vertex
         */
        int PopLowestVertex();
        /**
         * @brief rerturns the saturation degree corresponding to the GetHighestVertex()
         * 
         * @returns the saturation degree corresponding to the GetHighestVertex()
         */
        int GetHighestSatDegree() const;
        /**
         * @brief returns the vertex with the highest degree and highest saturation degree
         * 
         * @return the vertex
         */
        int GetHighestVertex() const;
        /**
         * @brief returns and remove the vertex with the highest degree and highest 
         *        saturation degree
         * 
         * @return the vertex
         */
        int PopHighestVertex();
        /**
         * @brief 
         * 
         */
        bool IsEmpty() const;
        /**
         * @brief returns a linked list of DSaturItems
         * 
         * @param sat_degree saturation degree corresponding to that list
         * @return DSaturItem* pointer to the beginning of the list
         */
        const DSaturItem* operator[](int sat_degree) const;

    protected:
        /**
         * @brief increases the saturation degree of the given vertex.
         * 
         * @param vertex vertex of which the saturation degree is increased
         * @returns true on success, else false
         */
        bool IncreaseSatDegree(int vertex, int increment=1);
        /**
         * @brief decreases the saturation degree of the given vertex.
         * 
         * @param vertex vertex of which the saturation degree is decreased
         * @returns true on success, else false
         */
        bool DecreaseSatDegree(int vertex, int decrement=-1);
        /**
         * @brief increases the degree of the given vertex.
         * 
         * @param vertex vertex of which the degree is increased
         * @returns true on success, else false
         */
        bool IncreaseDegree(int vertex, unsigned int increment=1);
        /**
         * @brief decreases the degree of the given vertex.
         * 
         * @param vertex vertex of which the degree is decreased
         * @param decrement a value > 0 which is the decrement of the degree
         * @returns true on success, else false
         */
        bool DecreaseDegree(int vertex, unsigned int decrement=1);


        /**
         * @brief maps each vertex to its DSaturItem
         */
        std::vector<DSaturItem*> _vertex_to_item;
        /**
         * @brief maps each vertex to the colors of its neighbours
         */
        std::vector<std::set<unsigned short>> _vertex_to_neighbour_colors;
        /**
         * @brief maps each saturation degree to its list of DSaturItems
         */
        std::vector<DSaturItem*> _sat_degree_to_list;
        DSaturItem* _last_item;
        unsigned int _last_degree;
};

#endif // DSATUR_COLOR_HPP

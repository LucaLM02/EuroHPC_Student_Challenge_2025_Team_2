#ifndef RECOLOR_HPP
#define RECOLOR_HPP

#include "graph.hpp"
#include <random>

/**
 *  @brief abstract functional class that wraps Recolor method; (partially) recolors a graph 
 *         to reduce the maximum color number
 *  @see the method description
 */
class RecolorStrategy {
    public:
        /**
         *  @brief tries to reduce the maximum color number of a certain coloring by 1
         * 
         *  @returns 0 if recoloring failed, otherwise how much the highest color was reduced
         */
        virtual unsigned int Recolor(Graph& graph) const = 0;
};

/**
 *  @brief functional class for the Recolor method; (partially) recolors a graph to reduce the maximum color number
 *  @see the method description
 */
class GreedySwapRecolorStrategy : public RecolorStrategy {
    public:
        /**
         *  @brief tries to reduce the maximum color number by 1 of a certain coloring by greedly 
         *         swapping colors
         * 
         *  @details For each vertex with the highest degree, reduces the color by swapping 
         *           it with a neighbour, which will be given another color. 
         * 
         * 
         *  @returns            0 if recoloring failed, otherwise how much the highest color was reduced
         */
        unsigned int Recolor(Graph& graph) const override;
};

/**
 * @brief contains all the data needed about a neigbhour of a vertex which 
 *        GreedySwapRecolorStrategy want to reduce the color
 */
class VertexRecolorData {
    public:
        const static unsigned short NOT_ASSIGNED = -1;

        VertexRecolorData();
        VertexRecolorData(int vertex,
                          std::vector<unsigned short>* coloring,
                          unsigned short max_color);

        /**
         * @brief initializes the vertex to which this data is assigned
         * 
         * @param current_color 
         */
        inline void InitVertex(int vertex, std::vector<int> neighbours) {
            _vertex = vertex;
            _neighbours = neighbours;
        }

        inline int GetVertex() {
            return _vertex;
        }

        /**
         * @brief gives to this class the global array in which the coloring of the graph
         *        is stored
         * 
         * @param colors 
         */
        inline void InitColoring(std::vector<unsigned short>* coloring, unsigned short max_color) {
            _coloring = coloring;
            _max_color = max_color;
        }

        inline unsigned short GetCurrentColor() {
            return (*_coloring)[_vertex];
        }

        /**
         * @brief reverts the color to the old one
         * 
         */
        inline void RevertColor() {
            if ( _old_color != NOT_ASSIGNED ) {
                (*_coloring)[_vertex] = _old_color;
                _old_color = NOT_ASSIGNED;
            }
        }

        /**
         * @brief whether to this vertex can be assigned a new color or not.
         * 
         * @return true 
         * @return false 
         */
        inline bool IsRecolorable() {  
            std::vector<bool> color_availability(_max_color+1, false);
            for ( int neighbour : _neighbours ) {
                color_availability[(*_coloring)[neighbour]] = false;

            }
            
            std::vector<int> available_colors;
            available_colors.reserve(_max_color+1);
            for (int i = 1; i < _max_color+1; i++) {
                if ( color_availability[i] ) {
                    available_colors.push_back(i);
                }
            }

            if ( color_availability.size() == 0 ) {
                return false;
            }
            return true;
        }

        /** 
         * @brief assigns a new color from the available ones.
         *        If the color is not the desired one, please do RevertColor() and
         *        then again thism method
        */
        inline bool Recolor() {

            if ( _neighbours.size() == 0 ) {
                return true;
            }

            // retrieves the available colors
            std::vector<bool> color_availability(_max_color, true);
            color_availability[this->GetCurrentColor()] = false;
            for ( int neighbour : _neighbours ) {
                color_availability[(*_coloring)[neighbour]] = false;

            }
            
            std::vector<int> available_colors;
            available_colors.reserve(_max_color);
            for (int i = 1; i < _max_color; i++) {
                if ( color_availability[i] ) {
                    available_colors.push_back(i);
                }
            }

            // not colorable
            if ( available_colors.size() == 0 ) {
                return false;
            } else if ( available_colors.size() == 1 ) {
                if ( _old_color == NOT_ASSIGNED ) {
                    _old_color = (*_coloring)[_vertex];
                }
                (*_coloring)[_vertex] = available_colors[0];
                return true;
            }


            if ( _old_color == NOT_ASSIGNED ) {
                _old_color = (*_coloring)[_vertex];
            }

            // randomly choosing a new color which is not forbidden
            std::uniform_int_distribution<int> u(0, available_colors.size()-1);
            (*_coloring)[_vertex] = available_colors[u(*_random_generator)];

            return true;
        }
        inline void AssignColor(unsigned short color);
        
    protected:
        int _vertex;
        /**
         * @brief previously assigned color, if any, otherwise NOT_ASSIGNED
         * 
         */
        unsigned short _old_color;
        /**
         * @brief available colors to swap with
         */
        unsigned short _max_color;
        std::vector<unsigned short>* _coloring;
        std::vector<int> _neighbours;
        std::unique_ptr<std::mt19937> _random_generator;
};

class SwapRecolorStructure {
    public:
        SwapRecolorStructure(Graph& graph, std::vector<unsigned short>& coloring, int threshold = 10);
        void FillWithData();

        bool Recolor();

    private:
        Graph& _graph;
        std::vector<VertexRecolorData> _vertex_to_data;
        std::set<int> _recolored_neighbours;
        std::vector<unsigned short>& _coloring;
        /**
         * @brief if the number of vertices to recolor (aka with maximum color) is
         *        greater or equal than threshold, do not color
         */
        int _threshold;
        bool _dont_color;

        bool RecolorBody(std::vector<int>& vertices);
};

#endif // RECOLOR_HPP

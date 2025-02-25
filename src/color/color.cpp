#include "color.hpp"

#include <cstring>
#include <iostream>
#include <cmath>

unsigned short GreedyFindColor(const Graph& graph,
                             const unsigned int vertex,
                             std::vector<unsigned short>& coloring, 
                             unsigned int current_max_k) {
    unsigned short neighbour_color;

    std::vector<int> neighbours;
    graph.GetNeighbours(vertex, neighbours);

    unsigned int max_colors = std::max(
                                static_cast<unsigned int>(neighbours.size()), 
                                current_max_k + 2);
    std::vector<bool> color_used(max_colors);   // a lower memory allocation should be performed
    //std::memset(&color_used[0], 0, color_used.size() * sizeof(color_used[0]));        
    std::fill(color_used.begin(), color_used.end(), false);

    // saves the color used by the neighbours, hence the non available ones
    for ( const int neighbour :  neighbours ) {
        neighbour_color = coloring[neighbour];
        if ( neighbour_color > 0 ) {
            color_used[neighbour_color - 1] = true;
        } else {
        }
    }


    // finds the first color available
    unsigned short k=0;
    while ( color_used[k] == true ) {
        k++;
    }
    for ( k = 0; color_used[k] == true; k++ ) {
    }

    return k+1;
}

void GreedyColorStrategy::Color(Graph& graph,
                                unsigned short& max_k) const
{ 
    // used for accessing neighbours colors
    std::vector<unsigned short> coloring(graph.GetHighestVertex() + 1, 0);

    graph.SortByDegree();

    unsigned short assigned_color;
    const std::vector<int>& original_vertices = graph.GetVertices();
    max_k = 0;

    // assigning a color to each vertex and storing it in its color class
    for ( int vertex : original_vertices ) {

        // probably can be done more efficiently but for now only readability is important
        assigned_color = GreedyFindColor(graph, vertex, coloring, max_k);
        if ( max_k < assigned_color ) {
            max_k = assigned_color;
        }

        coloring[vertex] = assigned_color;                      // saving color
    }

    graph.SetFullColoring(coloring);
}



#include "color.hpp"

#include <cstring>
#include <iostream>
#include <cmath>

unsigned int GreedyFindColor(const VertexSet& vertices,
                             const Edges& edges,
                             const unsigned int vertex_index,
                             std::vector<unsigned short>& coloring, 
                             unsigned int current_max_k,
                             const Constraints& constraints) {
    unsigned int neighbour_index;
    unsigned int neighbour_color;

    VertexSet neighbours;
    GetNeighbours(edges, vertex_index, neighbours);

    unsigned int max_colors = std::max(
                                static_cast<unsigned int>(neighbours.size()), 
                                current_max_k + 2);
    std::vector<bool> color_used(max_colors);   // a lower memory allocation should be performed
    //std::memset(&color_used[0], 0, color_used.size() * sizeof(color_used[0]));        
    std::fill(color_used.begin(), color_used.end(), false);

    // saves the color used by the neighbours, hence the non available ones
    for ( unsigned int neighbour_scroll; neighbour_scroll < neighbours.size(); neighbour_scroll++ ) {
        neighbour_index = neighbours[neighbour_scroll];

        neighbour_color = coloring[neighbour_index];
        if ( neighbour_color > 0 ) {
            color_used[neighbour_color - 1] = true;
        }
    }

    for ( unsigned int vertex_b : constraints.GetDiffConstraints(vertex_index)) {
        if ( coloring[vertex_b] != 0 ) {

        }
    }

    // finds the first color available
    int k=0;
    while ( color_used[k] == true ) {
        k++;
    }
    for ( k = 0; color_used[k] == true; k++ ) {
    }

    return k+1;
}

// TODO: config should be used
void GreedyColorStrategy::Color(VertexSet& vertices, const Edges& edges, 
                          const Constraints& constraints, 
                          std::vector<unsigned short>& coloring, 
                          unsigned short& max_k,
                          const unsigned int& expected_maximum_color) const
{ 
    // coloring algorithm explained in paper match2007 (An improved branch and bound algorithm for the maximum clique problem)
    // should be the best one
    // but actually here I do not have a k_min

    // maybe for the first recursion levels advanced algorithms such as evolutionary/genetic algorithms, tabu search, simulated annealing,
    // column generation could work

    // what about dsatur and recursive large first?

    std::vector<VertexSet> color_classes;
    color_classes.resize(expected_maximum_color);

    // used for accessing neighbours colors
    coloring.resize(vertices.size());
    std::fill(coloring.begin(), coloring.end(), 0);
    //std::memset(&local_coloring[0], 0, local_coloring.size() * sizeof(local_coloring[0]));        // more efficient

    max_k = 0;

    unsigned int vertex_index;
    unsigned short assigned_color;

    // assigning a color to each vertex and storing it in its color class
    for ( unsigned int vertex_scroll = 0; vertex_scroll < vertices.size(); vertex_scroll++ ) {
        vertex_index = vertices[vertex_scroll];

        // probably can be done more efficiently but for now only readability is important
        assigned_color = GreedyFindColor(vertices, edges, vertex_index, coloring, max_k, constraints);     
        if ( max_k < assigned_color ) {
            max_k = assigned_color;
        }

        coloring[vertex_index] = assigned_color;             // saving color
        color_classes[assigned_color-1].push_back(vertex_index);      // saving vertex into the color class
    }
        

    // transferring from color classes to vertices
    coloring.clear();
    vertices.clear();

    // remembering that colors starts from k=1...
    for ( unsigned short k=0; k<max_k; k++) {
        for (unsigned int vertex_scroll = 0; vertex_scroll < color_classes[k].size(); vertex_scroll++) {
            vertices.push_back(color_classes[k][vertex_scroll]);
            coloring.push_back(k+1);
        }
    }
}
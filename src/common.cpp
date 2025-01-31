#include "common.hpp"

void GetNeighbours(const Edges& edges, const unsigned int vertex_index, VertexSet& neighbours) { 
    neighbours.clear();
    neighbours.reserve(edges.size());    // it would be better density * num_vertices

    int counter = 0;
    for ( unsigned int i=0; i<edges.size(); i++ ) {
        if ( edges[vertex_index][i] ) {
            neighbours.push_back(i);
        }
    }
}
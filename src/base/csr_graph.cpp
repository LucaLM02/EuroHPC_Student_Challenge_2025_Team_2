#include "csr_graph.hpp"
#include <cstring>


// TODO: does the Dimacs contains both (v,w) and (w,v) ? no, only (v,w) -> ADAPT THE CODE
CSRGraph::CSRGraph(const Dimacs& dimacs_graph) {
    std::vector<Dimacs::Edge> unsorted_edges = dimacs_graph.edges;

    this->nVertices = dimacs_graph.numVertices;
    this->nEdges    = dimacs_graph.getNumEdges();

    this->offsets.reserve(this->nVertices);
    this->offsets[0] = 0;

    // offsets array from degrees
    for (int i=1; i<this->nVertices; i++) {
        this->offsets[i] = this->offsets[i-1] + dimacs_graph.degrees[i-1];
    }

    this->edges.reserve(this->nEdges);

    // given the offsets, building the edges vector
    std::vector<int> used_spaces(this->nVertices);
    std::fill(used_spaces.begin(), used_spaces.end(), 0);
    std::pair<int, int> tmp;
    for (int i=0; i<this->nEdges; i++) {
        tmp = dimacs_graph.edges[i];
        this->edges[used_spaces[tmp.first]] = tmp.second;
        used_spaces[tmp.first]++;
    }
}


CSRGraph::CSRGraph(const CSRGraph& other) {
    
}

/*
CSRGraph::CSRGraph(const CSRGraph& other, bool new_edge) {
    this->nVertices = other.nVertices;
    this->nEdges = other.nEdges;
    if ( new_edge ) 
        this->nEdges += 1;
    

}
*/

// TODO: ADD ALSO W-V EDGE
void CSRGraph::AddEdge(int v, int w) {
    
    if ( this->edges.size() == this->edges.max_size()) {
        this->edges.reserve(this->edges.size()*2);
    }

    // translates the elements after <v,w> edge
    this->edges.push_back(0);
    for (int i=this->nEdges-1; i>=this->offsets[v+1]; i--) {
        this->edges[i-1] = this->edges[i];
    }

    // adds <v,w> edge
    this->edges[this->offsets[v+1]-1] = w;

    // translated the offsets after <v,w> edge
    for (int i=v+1; i<this->nVertices; i++) {
        this->offsets[i]++;
    }

}

void CSRGraph::GetNeighbours(int vertex, std::vector<int> &result) const {
    result.clear();
}

void CSRGraph::GetNeighbours(int vertex, std::set<int> &result) const {
    result.clear();
}

size_t CSRGraph::GetNumVertices() const {
    return 0;
}

size_t CSRGraph::GetNumEdges() const {
    return 0;
}

void CSRGraph::MergeVertices(int v, int w) {
}

int CSRGraph::GetNeighboursIndex(int vertex) const {
    return 0;
}


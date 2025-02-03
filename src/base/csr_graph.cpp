#include "csr_graph.hpp"

#include <cstring>
#include <cmath>


// TODO: does the Dimacs contains both (v,w) and (w,v) ? no, only (v,w) -> ADAPT THE CODE
CSRGraph::CSRGraph(const Dimacs& dimacs_graph) {
    std::vector<Dimacs::Edge> unsorted_edges = dimacs_graph.edges;

    this->nVertices = dimacs_graph.numVertices;
    this->nEdges    = dimacs_graph.getNumEdges();

    this->offsets.resize(this->nVertices+1); // slot 0 is not used; vertices are counted from 1

    // offsets array from degrees
    for (int i=2; i<=this->nVertices; i++) {
        this->offsets[i] = this->offsets[i-1] + dimacs_graph.degrees[i-1];
    }

    this->edges.resize(this->nEdges);

    // given the offsets, building the edges vector
    std::vector<int> placed_edges_per_vertex(this->nVertices+1);
    std::fill(placed_edges_per_vertex.begin(), placed_edges_per_vertex.end(), 0);
    std::pair<int, int> edge;
    int edge_index;
    for (int i=0; i<this->nEdges; i++) {
        edge = dimacs_graph.edges[i];

        edge_index = this->offsets[edge.first] + placed_edges_per_vertex[edge.first];
        this->edges[edge_index] = edge.second;

        edge_index = this->offsets[edge.second] + placed_edges_per_vertex[edge.second];
        this->edges[edge_index] = edge.first;

        placed_edges_per_vertex[edge.first]++;
        placed_edges_per_vertex[edge.second]++;
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

void CSRGraph::AddEdge(int v, int w) {
    if ( this->edges.size() == this->edges.max_size()) {
        this->edges.reserve(this->edges.size()*2);
    }

    this->edges.push_back(0);
    this->edges.push_back(0);

    int max = std::max(v,w);
    const int min = std::min(v,w);

    nEdges++;

    // translates the elements after <max,min> vertex
    if ( v == w ) {
        // TODO
        _effective_nEdges++;    // only <v,w>=<v,v> added
    } else {

        /* offsets: |---A------B----|                       A = offset[min], B = offset[max]
         * edges:   |-----------(-------)-----[---]----|    () = neighbours of min, [] = neighbours of max
         *                                 
         * 1. The part ]----| needs to be translated by 2
         *      1.a translating edges
         *      1.b increasing by 2 the offsets
         * 2. Right after ] the element min is inserted    -> <max, min> edge
         * 3. The part )----] needs to be translated by 1
         *      3.a translating edges
         *      3.b increasing by 1 the offsets
         * 4. Right after ) the element max is inserted    -> <min, max> edge
         */

        int last_edge_of_max;
        if (  max < *_vertices.rbegin() ) {

            // 1.a
            for (int i = _effective_nEdges-1; i>=this->offsets[max+1]; i--) {
                this->edges[i+2] = this->edges[i];
            }

            // 2.a
            for (int i = max+1; i <= offsets.size(); i++) {
                offsets[i] += 2;
            }

            last_edge_of_max = this->offsets[max+1]-1;  // points to the first edge next to the last of max
        } else {
            // skipping 1.

            // 2 shifts so I write 3 position after (so +2)
            last_edge_of_max = _effective_nEdges + 2;       // points to the empty space next to the last edge of max
            max = *_vertices.rbegin();
        }
        
        // 2. adds <v,w> edge
        this->edges[last_edge_of_max] = min;

        // 3.a translates the elements between <min,max> and <max,min> vertex
        for (int i = last_edge_of_max-2; i>=this->offsets[min+1]; i--) {
            this->edges[i+1] = this->edges[i];
        }

        // 3.b translated the offsets after <v,w> edge
        for (int i=min+1; i<=max; i++) {
            this->offsets[i]++;
        }

        // 4. adds <v,w> edge
        this->edges[this->offsets[min+1]-1] = max;

        _effective_nEdges += 2;     // both <v,w> and <w,v> are added
    }
}

void CSRGraph::RemoveEdge(int v, int w) {
    // Opposite of what is done in AddEdge

    int x_index, y_index;
    int last_edge_of_max;

    int max = std::max(v,w);
    const int min = std::min(v,w);

    // translates the elements after <max,min> vertex
    if ( v == w ) {
        // TODO
        _effective_nEdges--;    // only <v,w>=<v,v> added
    } else {

        /* offsets: |---A------B----|                       A = offset[min], B = offset[max]
         * edges:   |-----------(---X---)-----[-Y-]----|    () = neighbours of min, [] = neighbours of max
         *                                 
         * 1. Searching where w is between ( and ): X
         * 2. Searching where v is between [ and ]: Y
         * 3. Translating elements from X to [ by Y
         *      3.a translating the edges by -1
         *      3.b shifting by -1 the offsets
         * 5. Translating elements from Y to |
         * 
         * note if v != w then min < max <= highest_vertex
         */


        // 1. finding X index
        x_index == -1;
        for (int i = this->offsets[min]; i < this->offsets[min+1]; i++) {
            if ( this->edges[i] == max ) {
                x_index = i;
                break;
            }
        }

        if ( x_index == -1 ) {
            return;
        }

        if ( max == *_vertices.rbegin() ) {
            last_edge_of_max = _effective_nEdges;
        } else {
            last_edge_of_max = this->offsets[max+1];
        }

        // 2. finding Y index
        for (int i = this->offsets[max]; i < last_edge_of_max; i++) {
            if ( this->edges[i] == min ) {
                y_index = i;
                break;
            }
        }

        // y_index != -1 for sure since x_index != -1

        // 3.a translates the elements between X and Y
        for (int i = x_index; i < y_index; i++) {
            this->edges[i] = this->edges[i+1];
        }

        _effective_nEdges--;    // one is removed

        // 3.b shifts by -1 the offsets between min and max
        for (int i=min+1; i<=max; i++) {
            this->offsets[i]--;
        }

        // 4.a translates the elements between Y and |
        for (int i = y_index+1; i < _effective_nEdges; i++) {
            this->edges[i-2] = this->edges[i];
        }

        _effective_nEdges--;     // also <max, min> is removed

        // 4.b shifts by -2 the offsets between max and the end
        for (int i = max+1; i < offsets.size(); i++) {
            this->offsets[i] -= 2;
        }

    }

}

void CSRGraph::AddVertex(int v) {
    if ( v <= 0 ) return;

    _vertices.insert(v);
    int next_v = *_vertices.upper_bound(v);
    
    offsets.push_back(0);
    // moving the offsets by 1 position
    for (int i=offsets.size()-1; i > next_v; i--) {
        offsets[i] = offsets[i-1];
    }

    offsets[v] = offsets[v-1];
}

void CSRGraph::RemoveVertex(int v) {
    if ( v < 0 || v > offsets.size() ) {
        return;
    }

    int delta;
    if ( v == *_vertices.rbegin() ) {
        delta = _effective_nEdges - this->offsets[v];
    } else {
        delta = this->offsets[v+1] - this->offsets[v];

        // removing edges of v
        std::memmove(&edges[this->offsets[v]], &edges[this->offsets[v+1]], _effective_nEdges - this->offsets[v+1]);

        // removing also all the other vertices
    }

    // decreasing offsets
    for ( int i=v+1; i < offsets.size(); i++ ) {
        offsets[i] -= delta;
    }

    _effective_nEdges -= delta*2;

    _vertices.erase(v);
}

// TODO: finish this
void CSRGraph::RemoveVertexWithRenaming(int v) {
    this->RemoveVertex(v);
}

void CSRGraph::GetNeighbours(int vertex, std::vector<int> &result) const {
    result.clear();

    int end_edge_index;
    if ( vertex == *_vertices.rbegin() ) {
        end_edge_index = _effective_nEdges;
    } else {
        end_edge_index = this->offsets[vertex+1];
    }

    result.reserve(end_edge_index - this->offsets[vertex]);

    for ( int i = this->offsets[vertex]; i < end_edge_index; i++ ) {
        result.push_back(edges[i]);
    }
}

void CSRGraph::GetNeighbours(int vertex, std::set<int> &result) const {
    result.clear();

    int end_edge_index;
    if ( vertex == *_vertices.rbegin() ) {
        end_edge_index = _effective_nEdges;
    } else {
        end_edge_index = this->offsets[vertex+1];
    }

    for ( int i = this->offsets[vertex]; i < end_edge_index; i++ ) {
        result.insert(edges[i]);
    }
}

size_t CSRGraph::GetNumVertices() const {
    return _vertices.size();
}

size_t CSRGraph::GetNumEdges() const {
    return 
}

void CSRGraph::MergeVertices(int v, int w) {
}

int CSRGraph::GetNeighboursIndex(int vertex) const {
    return 0;
}


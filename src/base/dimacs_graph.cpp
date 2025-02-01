#include "dimacs_graph.hpp"

DimacsGraph::DimacsGraph(Dimacs& dimacs) 
: _dimacs{dimacs}
{
    // NOTE: assuming _dimacs has continuous vertices from i=0 to i=numVertices-1
    for ( int i=0; i<_dimacs.numVertices; i++) {
        _vertices.insert(i);
    }
}

void DimacsGraph::AddEdge(int v, int w) {
    _dimacs.edges.push_back(std::make_pair(v, w));
    _dimacs.degrees[v]++;
    _dimacs.degrees[w]++;
}

void DimacsGraph::RemoveEdge(int v, int w) {
    std::pair<int, int> edge;
    for (int i = 0; i < _dimacs.edges.size(); i++) {
        edge = _dimacs.edges[i];

        if ( (edge.first == v && edge.second == w) ||
             (edge.first == w && edge.second == v)) {
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
                return;
        }
    }
}

void DimacsGraph::AddVertex(int v) {
    _vertices.insert(v);
}

void DimacsGraph::RemoveVertex(int v) {
    _vertices.erase(v);
}

void DimacsGraph::GetNeighbours(int vertex, std::vector<int> &result) const {
    result.clear();
    for ( const std::pair<int, int> &edge : _dimacs.edges ) {
        if ( edge.first == vertex )
            result.push_back(edge.second);
        else if (edge.second == vertex )
            result.push_back(edge.first);
    }
}

void DimacsGraph::GetNeighbours(int vertex, std::set<int> &result) const {
    result.clear();
    for ( const std::pair<int, int> &edge : _dimacs.edges ) {
        if ( edge.first == vertex )
            result.insert(edge.second);
        else if (edge.second == vertex )
            result.insert(edge.first);
    }
}

bool DimacsGraph::HasEdge(int v, int w) const {
    for ( std::pair<int, int> edge : _dimacs.edges ) {
        if ( (edge.first == v && edge.second == w) ||
             (edge.first == w && edge.second == v)) {
                return true;
             }
    }
    return false;
}

void DimacsGraph::GetEdges(int vertex, std::set<std::pair<int,int>> &result) const {

}

void DimacsGraph::GetEdges(int vertex, std::vector<std::pair<int, int>> &result) const {

}

size_t DimacsGraph::GetNumVertices() const {
    return _dimacs.numVertices;
}

size_t DimacsGraph::GetNumEdges() const {
    return _dimacs.getNumEdges();
}

bool DimacsGraph::MergeVertices(int v, int w) {

    std::set<int> v_neighbours;
    this->GetNeighbours(v, v_neighbours);

    std::pair<int, int> edge;
    for ( int i=0; i < _dimacs.edges.size(); i++ ) {
        edge = _dimacs.edges[i];
        if (edge.first == w) {
            if (!v_neighbours.contains(edge.second)) 
                // adding vertex as neighbour of v
                _dimacs.edges[i].first = v;
            else 
                // duplicate, erasing it
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
        } else if (edge.second == w) {
            if (!v_neighbours.contains(edge.first))
                // adding vertex as neighbour of v
                _dimacs.edges[i].second = v;
            else
                // duplicate, erasing it
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
        } 
    }

    // "removing" vertex v
    _dimacs.numVertices--;      // however the vertices are not contiguous necessarly!
}
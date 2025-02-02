#include "dimacs_graph.hpp"

#include <iostream>

DimacsGraph::DimacsGraph(Dimacs& dimacs) 
: _dimacs{dimacs}
{
    // NOTE: assuming _dimacs has continuous vertices from i=1 to i=numVertices
    for ( int i=1; i<=_dimacs.numVertices; i++) {
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
    bool has_removed=false;
    for (int i = 0; i < _dimacs.edges.size(); i++) {
        edge = _dimacs.edges[i];
        if ( has_removed ) {
            i--;
            has_removed = false;
        }

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
    _dimacs.numVertices--;      // however the vertices are not contiguous necessarly!

    std::pair<int, int> edge;
    bool has_removed = false;
    for (int i = 0; i < _dimacs.edges.size(); i++) {
        if ( has_removed ) {
            i--;
            has_removed = false;
        }
        edge = _dimacs.edges[i];
        if ( edge.first == v || edge.second == v ) {
            _dimacs.edges.erase(_dimacs.edges.begin() + i);
            has_removed = true;
        }
    }
}

// TODO: finish this
void DimacsGraph::RemoveVertexWithRenaming(int v) {
    this->RemoveVertex(v);
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

void DimacsGraph::GetVertices(std::set<int> &result) const {
    result = this->_vertices;
}
void DimacsGraph::GetVertices(std::vector<int> &result) const {
    std::copy(_vertices.begin(), _vertices.end(), std::back_inserter(result));
}

size_t DimacsGraph::GetNumVertices() const {
    return _dimacs.numVertices;
}

size_t DimacsGraph::GetNumEdges() const {
    return _dimacs.getNumEdges();
}

void DimacsGraph::MergeVertices(int v, int w) {

    std::set<int> v_neighbours;
    this->GetNeighbours(v, v_neighbours);

    std::cout << "Neighbours of " << v << ":";
    std::cout << std::endl;

    std::pair<int, int> edge;
    bool has_removed = false;
    for ( int i=0; i < _dimacs.edges.size(); i++ ) {
        if ( has_removed ) {
            i--;
            has_removed = false;
        }
        edge = _dimacs.edges[i];
        if (edge.first == w) {
            if (!v_neighbours.contains(edge.second)) {
                // adding vertex as neighbour of v
                _dimacs.edges[i].first = v;
            } else {
                // duplicate, erasing it
                has_removed = true;
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
            }
        } else if (edge.second == w) {
            if (!v_neighbours.contains(edge.first)) {
                // adding vertex as neighbour of v
                _dimacs.edges[i].second = v;
            } else {
                // duplicate, erasing it
                has_removed = true;
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
            }
        } 
    }

    // "removing" vertex v
    has_removed = false;
    for (int i=0; i<_dimacs.edges.size(); i++) {
        if ( has_removed ) {
            i--;
            has_removed = false;
        }

        if ( _dimacs.edges.size() == w ) 
            _dimacs.edges.erase(_dimacs.edges.begin()+i);
    }
    this->RemoveVertex(w);

}
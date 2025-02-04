#include "csr_graph.hpp"

#include <iostream>
#include <cstring>
#include <cmath>

CSRGraph* CSRGraph::LoadFromDimacs(const std::string& file_name) {
    Dimacs dimacs;
    dimacs.load(file_name.c_str());
    CSRGraph* graph = new CSRGraph(dimacs);
    return graph;
}

void CSRGraph::AddEdge(int v, int w) {
    _edges[v].push_back(w);
    _edges[w].push_back(v);
    _nEdges++;
    _cache_degrees_invalidated = true;
}

void CSRGraph::RemoveEdge(int v, int w) {
    // efficiently removes the edges by taking the edge to remove, swapping it with the 
    // last element and then removing it
    // complexity is still O(n)
    auto it = std::find(_edges[v].begin(), _edges[v].end(), w);
    if (it != _edges[v].end()) {
        std::swap(*it, _edges[v].back()); 
        _edges[v].pop_back(); 
        _cache_degrees_invalidated = true;
        _nEdges--;
    }

    it = std::find(_edges[w].begin(), _edges[w].end(), v);
    if (it != _edges[w].end()) {
        std::swap(*it, _edges[w].back()); 
        _edges[w].pop_back(); 
        _cache_degrees_invalidated = true;
    }


}

int CSRGraph::AddVertex() {
    int v = _vertices.size();

    _vertices.push_back(v);
    _edges.emplace_back(0);

    return v;
}

void CSRGraph::RemoveVertex(int v) {
    // removing the vertex
    _vertices.erase(std::find(_vertices.begin(), _vertices.end(), v));

    // removing the vertices
    _nEdges -= _edges[v].size();
    _edges[v].clear();
    for ( int vertex : _vertices ) {
        auto it = std::find(_edges[vertex].begin(), _edges[vertex].end(), v);
        if (it != _edges[vertex].end()) {
            std::swap(*it, _edges[vertex].back()); 
            _edges[vertex].pop_back(); 
            _cache_degrees_invalidated = true;
        }
    }


}

// TODO: finish this
void CSRGraph::RemoveVertexWithRenaming(int v) {
    this->RemoveVertex(v);
}


void CSRGraph::MergeVertices(int v, int w) {
    // O(2*#neighbours*log(#neighbours))
    std::sort(_edges[v].begin(), _edges[v].end());
    std::sort(_edges[w].begin(), _edges[w].end());


    std::vector<int> deleted_edges;
    std::vector<int> modified_edges;
    deleted_edges.reserve(_edges[w].size());
    modified_edges.reserve(_edges[w].size());

    _edges[v].reserve(_edges[v].size() + _edges[w].size());

    // merging neighbours of w into v, avoiding duplicates
    // O(#neighbours*log(#neighbours))
    for ( const int nw : _edges[w] ) {
        if ( !std::binary_search(_edges[v].begin(), _edges[v].end(), nw) ) {
            _edges[v].push_back(nw);
            modified_edges.push_back(nw);
        } else {
            deleted_edges.push_back(nw);
        }
    }

    _edges[w].clear();
    _vertices.erase(std::find(_vertices.begin(), _vertices.end(), w));

    // deleting `w` from the neighbour lists of common neighbours between `v` and `w`
    for ( const int deleted_edge : deleted_edges ) {

        // O(#neighbours) but I rarely iterate over the full vector
        auto it = std::find(_edges[deleted_edge].begin(), _edges[deleted_edge].end(), w);
        if (it != _edges[deleted_edge].end()) {
            std::swap(*it, _edges[deleted_edge].back()); 
            _edges[deleted_edge].pop_back(); 
            _cache_degrees_invalidated = true;
            _nEdges--;
        }
    }

    // modifying `w` into `v` in the neighbour lists of all the other neighbours of `w`
    for ( const int modified_edge : modified_edges ) {

        // O(#neighbours) but I rarely iterate over the full vector
        for (int i = 0; i < _edges[modified_edge].size(); i++) {
            if ( _edges[modified_edge][i] == w ) {
                _edges[modified_edge][i] = v;
                break;
            }
        }
    }

}

void CSRGraph::GetNeighbours(int vertex, std::vector<int> &result) const {
    result = _edges[vertex];
}

void CSRGraph::GetNeighbours(int vertex, std::set<int> &result) const {
    for ( int w : _edges[vertex] ) {
        result.insert(w);
    }
    
}

bool CSRGraph::HasEdge(int v, int w) const {
    // with this check I search through the shorter vector
    if ( _edges[v].size() > _edges[w].size() ) {
        return ( std::find(_edges[w].begin(), _edges[w].end(), v) != _edges[v].end());
    } else {
        return ( std::find(_edges[v].begin(), _edges[v].end(), w) != _edges[v].end());
    }
}

void CSRGraph::GetUnorderedVertices(std::set<int> &result) const {
    for ( int vertex : _vertices ) {
        result.insert(vertex);
    }
}

const std::vector<int>& CSRGraph::GetVertices() const {
    return _vertices;
}

void CSRGraph::SetVertices(std::vector<int>& vertices) {
    _vertices = vertices;
    _cache_degrees_invalidated = true;
}

size_t CSRGraph::GetNumVertices() const {
    return _vertices.size();
}

size_t CSRGraph::GetNumEdges() const {
    return _nEdges;
}

unsigned int CSRGraph::GetDegree(int vertex) const {
    return _edges[vertex].size();
}

const std::vector<int>& CSRGraph::GetDegrees() const {
    if ( _cache_degrees_invalidated ) {
        _ComputeCacheDegrees();
        _cache_degrees_invalidated = false;
    }
    return _cache_degrees;
}

unsigned int CSRGraph::GetMaxDegree() const {
    if ( _cache_degrees_invalidated ) {
        _ComputeCacheDegrees();
        _cache_degrees_invalidated = false;
    }

    return *std::max_element(_cache_degrees.begin(), _cache_degrees.end());
}

int CSRGraph::GetVertexWithMaxDegree() const {
    if ( _cache_degrees_invalidated ) {
        _ComputeCacheDegrees();
        _cache_degrees_invalidated = false;
    }

    int max_index = std::distance(_cache_degrees.begin(), 
                                  std::max_element(
                                    _cache_degrees.begin(), 
                                    _cache_degrees.end()));
    return _vertices[max_index];
}

// ------------------------ PROTECTED --------------------------
CSRGraph::CSRGraph(const Dimacs& dimacs_graph) 
: _vertices(static_cast<int>(dimacs_graph.numVertices)), 
  _nEdges{dimacs_graph.getNumEdges()},
  _edges(dimacs_graph.numVertices + 1u)
{

    int size = _vertices.size();
    for ( int vertex = 1; vertex <= size; vertex++ ) {
        _vertices[vertex-1] = vertex;
        _edges[vertex].reserve(dimacs_graph.degrees[vertex]);
    }

    for ( const std::pair<int, int>& edge : dimacs_graph.edges ) {
        _edges[edge.first].push_back(edge.second);
        _edges[edge.second].push_back(edge.first);
    }

}

void CSRGraph::_ComputeCacheDegrees() const {
    _cache_degrees.clear();
    _cache_degrees.resize(_vertices.size());
    for (int i = 0; i < _vertices.size(); i++) {
        _cache_degrees[i] = _edges[_vertices[i]].size();
    }
}
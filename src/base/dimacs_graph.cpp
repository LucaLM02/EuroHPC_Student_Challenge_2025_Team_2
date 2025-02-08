#include "dimacs_graph.hpp"

#include <iostream>
#include <algorithm>

DimacsGraph* DimacsGraph::LoadFromDimacs(const std::string& file_name) {
    return new DimacsGraph(file_name);
}

void DimacsGraph::AddEdge(int v, int w) {
    if ( v >= _dimacs.degrees.size() || w >= _dimacs.degrees.size() ) {
        return;
    }

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
                _dimacs.degrees[v]--;
                _dimacs.degrees[w]--;
                return;
        }
    }
    
}

int DimacsGraph::AddVertex() {
    int v = _vertices.size();
    _vertices.push_back(v);

    _dimacs.degrees.push_back(0);
    return v;
}

void DimacsGraph::RemoveVertex(int v) {
    _RemoveOnlyVertex(v);

    std::pair<int, int> edge;
    bool has_removed = false;
    for (int i = 0; i < _dimacs.edges.size(); i++) {
        if ( has_removed ) {
            i--;
            has_removed = false;
        }
        edge = _dimacs.edges[i];
        if ( edge.first == v ) { 
            _dimacs.degrees[edge.second]--;
            _dimacs.edges.erase(_dimacs.edges.begin() + i);
            has_removed = true;
        } 
        else if ( edge.second == v ) {
            _dimacs.degrees[edge.first]--;
            _dimacs.edges.erase(_dimacs.edges.begin() + i);
            has_removed = true;
        }
    }
    _RemoveDegree(v);
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

void DimacsGraph::GetUnorderedVertices(std::set<int> &result) const {
    result = std::set<int>(_vertices.begin(), _vertices.end());

}

const std::vector<int>& DimacsGraph::GetVertices() const {
    return _vertices;
}

int DimacsGraph::GetVertexByIndex(int index) const
{
    return _vertices[index];
}

int DimacsGraph::GetHighestVertex() const
{
    int max_vertex = 0;
    for ( int vertex : _vertices ) {
        if ( vertex > max_vertex ) {
            max_vertex = vertex;
        }
    }
    return max_vertex;
}

const std::set<int> &DimacsGraph::GetDeletedVertices() const
{
    return _deleted_vertices;
}

void DimacsGraph::SetVertices(std::vector<int> &vertices) {
    this->_vertices = vertices;
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
    unsigned int delta_degree = 0;
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
                delta_degree++;
            } else {
                // duplicate, erasing it
                has_removed = true;
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
                _dimacs.degrees[_dimacs.edges[i].second]--;     // vertex `second` looses 1 degree
            }
        } else if (edge.second == w) {
            if (!v_neighbours.contains(edge.first)) {
                // adding vertex as neighbour of v
                _dimacs.edges[i].second = v;
                delta_degree++;
            } else {
                // duplicate, erasing it
                has_removed = true;
                _dimacs.edges.erase(_dimacs.edges.begin()+i);
                _dimacs.degrees[_dimacs.edges[i].first]--;      // vertex `second` looses 1 degree
            }
        } 
    }

    _dimacs.degrees[v] += delta_degree;

    // "removing" vertex w
    this->_RemoveOnlyVertex(w);
}


unsigned int DimacsGraph::GetDegree(int vertex) const {
    return _dimacs.degrees[vertex];
}
std::vector<int> DimacsGraph::GetDegrees() const {
    // horrible and stupid solution, it's temporary
    _tmp_degrees.clear();
    _tmp_degrees.reserve(_vertices.size());
    for ( int vertex : _vertices ) {
        _tmp_degrees.push_back(_dimacs.degrees[vertex]);
    }

    return _tmp_degrees;
}

void DimacsGraph::GetDegrees(std::vector<int> &result) const
{
    result.clear();
    result.reserve(_vertices.size());
    for ( int vertex : _vertices ) {
        result.push_back(_dimacs.degrees[vertex]);
    }

}

unsigned int DimacsGraph::GetMaxDegree() const {
    return *std::max_element(_dimacs.degrees.begin(), _dimacs.degrees.end());
}

int DimacsGraph::GetVertexWithMaxDegree() const {
    return std::distance(
                _dimacs.degrees.begin(), 
                std::max_element(_dimacs.degrees.begin(), _dimacs.degrees.end())
           );
}

std::unique_ptr<Graph> DimacsGraph::Clone() const
{
    std::unique_ptr<DimacsGraph> clone = std::make_unique<DimacsGraph>(*this);

    return std::move(clone);
}

// ---------------------------- PROTECTED --------------------------------
DimacsGraph::DimacsGraph(const std::string& file_name) {
    _dimacs.load(file_name.c_str());
    // NOTE: assuming _dimacs has continuous vertices from i=1 to i=numVertices
    for ( int i=1; i<=_dimacs.numVertices; i++) {
        _vertices.push_back(i);
    }

}

// ----------------------------- PRIVATE ----------------------------------
void DimacsGraph::_RemoveOnlyVertex(int v) {
    _deleted_vertices.insert(v);
    _vertices.erase(std::find(_vertices.begin(), _vertices.end(), v));
    _dimacs.numVertices--;      // however the vertices are not contiguous necessarly!
}

void DimacsGraph::_RemoveDegree(int vertex) {
    _dimacs.degrees[vertex] = 0;
}

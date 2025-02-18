#include "csr_graph.hpp"

#include <cmath>
#include <cstring>

CSRGraph* CSRGraph::LoadFromDimacs(const std::string& file_name) {
	Dimacs dimacs;
	dimacs.load(file_name.c_str());
	CSRGraph* graph = new CSRGraph(dimacs);
	return graph;
}

CSRGraph::CSRGraph()
    : _vertices(0),
      _degrees(1),
      _coloring(1),
      _edges(1),
      _nEdges(0),
      _max_vertex(0) {}

	  bool CSRGraph::isEqual(const Graph &ot) const
	  {
		  CSRGraph other = dynamic_cast<const CSRGraph&>(ot);
		  if ( this->_vertices.size() != other._vertices.size() ) {
			  std::cout << "Different number of vertices: " << _vertices.size() << " vs "
						<< other._vertices.size() << std::endl;
	  
			  for ( int vertex : _vertices ) {
				  if ( std::find(other._vertices.begin(), other._vertices.end(), vertex) 
						  == other._vertices.end() ) {
					  std::cout << "Vertex: " << vertex << " is not contained in other graph" << std::endl;
					  continue;
				  }
				  if ( _edges[vertex].size() != other._edges[vertex].size() ) {
					  std::cout << "Edge lists do not have same size: " << _edges[vertex].size() << " vs "
								<< other._edges[vertex].size() << std::endl;
				  }
				  for ( int edge : _edges[vertex] ) {
					  if ( std::find(other._edges[vertex].begin(), other._edges[vertex].end(), edge) 
							  == other._edges[vertex].end() ) {
						  std::cout << "Edge : " << edge << " is not contained in other graph" << std::endl;
					  }
				  }
			  }
	  
		  }
		  for ( int vertex : _vertices ) {
			if ( std::find(other._vertices.begin(), other._vertices.end(), vertex) 
					== other._vertices.end() ) {
				std::cout << "Vertex: " << vertex << " is not contained in other graph" << std::endl;
				continue;
			}
			if ( _edges[vertex].size() != other._edges[vertex].size() ) {
				std::cout << "Edge lists do not have same size: " << _edges[vertex].size() << " vs "
						  << other._edges[vertex].size() << std::endl;
			}
			for ( int edge : _edges[vertex] ) {
				if ( std::find(other._edges[vertex].begin(), other._edges[vertex].end(), edge) 
						== other._edges[vertex].end() ) {
					std::cout << "Edge : " << edge << " is not contained in other graph" << std::endl;
				}
			}
			if(_degrees.size() != other._degrees.size()){
				std::cout << "Degrees size is different" << _degrees.size() << " " << other._degrees.size() << std::endl;
			}
			for(int i = 0; i < _degrees.size(); i++){
				if(_degrees[i] != other._degrees[i]){
					std::cout << "Degrees are different" << std::endl;
				}
			}
		}
		/*
		if(_coloring.size() != other._coloring.size()){
			std::cout << "Coloring size is different " << _coloring.size() << " " << other._coloring.size() << std::endl;
		}
		for(int i = 0; i < _coloring.size(); i++){
			if(_coloring[i] != other._coloring[i]){
				std::cout << "Coloring is different" << _coloring[i] << " " << other._coloring[i] <<  std::endl;
			}
		}
			*/
		  return true;
	  } 

std::string CSRGraph::Serialize() const {
	std::ostringstream oss;
	oss << _vertices.size() << " " << _nEdges << " " << _max_vertex << "\n";
	for (int vertex : _vertices) {
		oss << vertex << " ";
	}
	oss << "\n";
	for (int degree : _degrees) {
		oss << degree << " ";
	}
	oss << "\n";
	for (const auto& edges : _edges) {
		oss << edges.size() << " ";
		for (int edge : edges) {
			oss << edge << " ";
		}
	}
	oss << "\n";
	for (unsigned short color : _coloring) {
		oss << color << " ";
	}
	oss << "\n";
	return oss.str();
}

void CSRGraph::Deserialize(const std::string& data) {
	std::istringstream iss(data);
	size_t numVertices, numEdges;
	iss >> numVertices >> numEdges >> _max_vertex;

	_vertices.resize(numVertices);
	for (size_t i = 0; i < numVertices; ++i) {
		iss >> _vertices[i];
	}

	_degrees.resize(numVertices + 1);
	for (size_t i = 0; i <= numVertices; ++i) {
		iss >> _degrees[i];
	}

	_edges.resize(numVertices + 1);
	for (size_t i = 0; i <= numVertices; ++i) {
		int edgeSize;
		iss >> edgeSize;
		for(int j = 0; j < edgeSize; ++j) {
			int edge;
			iss >> edge;
			_edges[i].push_back(edge);
		}
	}

	_coloring.resize(numVertices + 1);
	for (size_t i = 0; i <= numVertices; ++i) {
		iss >> _coloring[i];
	}

	_nEdges = numEdges;
}

void CSRGraph::AddEdge(int v, int w) {
	_edges[v].push_back(w);
	_edges[w].push_back(v);
	_nEdges++;

	_degrees[v]++;
	_degrees[w]++;
}

void CSRGraph::RemoveEdge(int v, int w) {
	// efficiently removes the edges by taking the edge to remove, swapping
	// it with the last element and then removing it complexity is still
	// O(n)
	auto it = std::find(_edges[v].begin(), _edges[v].end(), w);
	if (it != _edges[v].end()) {
		std::swap(*it, _edges[v].back());
		_edges[v].pop_back();
		_degrees[v]--;
		_nEdges--;
	}

	it = std::find(_edges[w].begin(), _edges[w].end(), v);
	if (it != _edges[w].end()) {
		std::swap(*it, _edges[w].back());
		_edges[w].pop_back();
		_degrees[w]--;
	}
}

int CSRGraph::AddVertex() {
	int v = _max_vertex + 1;
	_max_vertex++;

	_vertices.push_back(v);
	_degrees.emplace_back(0);
	_coloring.emplace_back(0);
	_edges.emplace_back(0);

	return v;
}

void CSRGraph::RemoveVertex(int v) {
	// removing the vertex
	_vertices.erase(std::find(_vertices.begin(), _vertices.end(), v));
	_removed_vertices.insert(v);

	// removing the vertices
	_nEdges -= _edges[v].size();
	_edges[v].clear();
	for (int vertex : _vertices) {
		auto it =
		    std::find(_edges[vertex].begin(), _edges[vertex].end(), v);
		if (it != _edges[vertex].end()) {
			std::swap(*it, _edges[vertex].back());
			_edges[vertex].pop_back();
			_degrees[vertex]--;
		}
	}

	_degrees[v] = 0;
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
        if ( !std::binary_search(_edges[v].begin(), _edges[v].end(), nw) && nw != v ) {
            _edges[v].push_back(nw);
            modified_edges.push_back(nw);
            _degrees[v]++;
        } else {
            deleted_edges.push_back(nw);
        }
    }

    _edges[w].clear();
    _vertices.erase(std::find(_vertices.begin(), _vertices.end(), w));
    _removed_vertices.insert(w);
    _degrees[w] = 0;

    // deleting `w` from the neighbour lists of common neighbours between `v` and `w`
    for ( const int deleted_edge : deleted_edges ) {

        // O(#neighbours) but I rarely iterate over the full vector
        auto it = std::find(_edges[deleted_edge].begin(), _edges[deleted_edge].end(), w);
        if (it != _edges[deleted_edge].end()) {
            std::swap(*it, _edges[deleted_edge].back()); 
            _edges[deleted_edge].pop_back(); 
            _degrees[deleted_edge]--;
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

void CSRGraph::SetColoring(const std::vector<unsigned short>& colors)
{
    _coloring;
    for (int i = 0; i < _vertices.size(); i++ ) {
        _coloring[_vertices[i]] = colors[i];
    }
}

void CSRGraph::SetColoring(int vertex, unsigned short color)
{
    _coloring[vertex] = color;
}

void CSRGraph::SetFullColoring(const std::vector<unsigned short> &colors)
{
    _coloring = colors;
}

void CSRGraph::ClearColoring()
{
    for (int vertex : _vertices ) {
        _coloring[vertex] = 0;
    }
}

void CSRGraph::SortByDegree(bool ascending)
{
    auto ascendingCompare = 
    [&](int v, int w) -> bool {
        return _degrees[v] < _degrees[w];
    };
    auto descendingCompare = 
    [&](int v, int w) -> bool {
        return _degrees[v] > _degrees[w];
    };

    if ( ascending ) {
        std::sort(_vertices.begin(), _vertices.end(), ascendingCompare);
    } else {
        std::sort(_vertices.begin(), _vertices.end(), descendingCompare);
    }
}

void CSRGraph::SortByExDegree(bool ascending)
{
    std::vector<int> ex_degrees(_degrees.size());
    std::vector<int> neighbours;
    for ( int vertex : _vertices ) {
        ex_degrees[vertex] = GetExDegree(vertex);
    }


    auto ascendingCompare = 
    [&](int v, int w) -> bool {
        return ex_degrees[v] < ex_degrees[w];
    };
    auto descendingCompare = 
    [&](int v, int w) -> bool {
        return ex_degrees[v] > ex_degrees[w];
    };

    if ( ascending ) {
        std::sort(_vertices.begin(), _vertices.end(), ascendingCompare);
    } else {
        std::sort(_vertices.begin(), _vertices.end(), descendingCompare);
    }

}

void CSRGraph::SortByColor(bool ascending)
{
    static auto ascendingCompare = 
    [&](int v, int w) -> bool {
        return _coloring[v] < _coloring[w];
    };
    static auto descendingCompare = 
    [&](int v, int w) -> bool {
        return _coloring[v] > _coloring[w];
    };

    if ( ascending ) {
        std::sort(_vertices.begin(), _vertices.end(), ascendingCompare);
    } else {
        std::sort(_vertices.begin(), _vertices.end(), descendingCompare);
    }
}

void CSRGraph::GetNeighbours(int vertex, std::vector<int> &result) const {
    result.assign(_edges[vertex].begin(), _edges[vertex].end());
}

void CSRGraph::GetNeighbours(int vertex, std::set<int> &result) const {
    result.clear();
    for ( int w : _edges[vertex] ) {
        result.insert(w);
    }
    
}

bool CSRGraph::HasEdge(int v, int w) const {
	// with this check I search through the shorter vector
	if (_edges[v].size() > _edges[w].size()) {
		return (std::find(_edges[w].begin(), _edges[w].end(), v) !=
			_edges[w].end());
	} else {
		return (std::find(_edges[v].begin(), _edges[v].end(), w) !=
			_edges[v].end());
	}
}

void CSRGraph::GetUnorderedVertices(std::set<int>& result) const {
	for (int vertex : _vertices) {
		result.insert(vertex);
	}
}

const std::vector<int>& CSRGraph::GetVertices() const { return _vertices; }

int CSRGraph::GetVertexByIndex(int index) const { return _vertices[index]; }

int CSRGraph::GetHighestVertex() const {
	int max_vertex = 0;
	for (int vertex : _vertices) {
		if (vertex > max_vertex) {
			max_vertex = vertex;
		}
	}
	return max_vertex;
}

const std::set<int>& CSRGraph::GetDeletedVertices() const {
	return _removed_vertices;
}

void CSRGraph::SetVertices(std::vector<int>& vertices) { _vertices = vertices; }

size_t CSRGraph::GetNumVertices() const { return _vertices.size(); }

size_t CSRGraph::GetNumEdges() const { return _nEdges; }

unsigned int CSRGraph::GetDegree(int vertex) const {
    return _degrees[vertex];
}

std::vector<int> CSRGraph::GetDegrees() const {
	std::vector<int> degrees(_vertices.size());
	this->GetDegrees(degrees);
	return degrees;
}

std::vector<int> CSRGraph::GetFullDegrees() const { return _degrees; }

void CSRGraph::GetFullDegrees(std::vector<int>& result) const {
	result = _degrees;
}

void CSRGraph::GetDegrees(std::vector<int>& result) const {
	result.clear();
	result.reserve(_degrees.size());
	for (int i = 0; i < _vertices.size(); i++) {
		result.push_back(_degrees[_vertices[i]]);
	}
}

unsigned int CSRGraph::GetMaxDegree() const {
	return *std::max_element(_degrees.begin(), _degrees.end());
}

int CSRGraph::GetVertexWithMaxDegree() const {
	int max_index =
	    std::distance(_degrees.begin(),
			  std::max_element(_degrees.begin(), _degrees.end()));
	return max_index;
}

int CSRGraph::GetExDegree(int vertex) const {
	int ex_degree = 0;
	for (int neighbour : _edges[vertex]) {
		ex_degree += _degrees[vertex];
	}

	return ex_degree;
}

std::vector<int> CSRGraph::GetMergedVertices(int vertex) const { return {}; }

std::vector<unsigned short> CSRGraph::GetColoring() const {
	std::vector<unsigned short> colors(_vertices.size());

	for (int i = 0; i < _vertices.size(); i++) {
		colors[i] = _coloring[_vertices[i]];
	}

	return colors;
}

std::vector<unsigned short> CSRGraph::GetFullColoring() const {
	return _coloring;
}

unsigned short CSRGraph::GetColor(int vertex) const {
	return _coloring[vertex];
}

std::unique_ptr<Graph> CSRGraph::Clone() const {
	std::unique_ptr<CSRGraph> graph = std::make_unique<CSRGraph>(*this);

	return std::move(graph);
}

// ------------------------ PROTECTED --------------------------
CSRGraph::CSRGraph(const Dimacs& dimacs_graph) 
: _vertices(static_cast<int>(dimacs_graph.numVertices)), 
  _nEdges{dimacs_graph.getNumEdges()},
  _edges(dimacs_graph.numVertices + 1u),
  _coloring(dimacs_graph.numVertices + 1u),
  _max_vertex(dimacs_graph.numVertices)
{

    int size = _vertices.size();
    for ( int vertex = 1; vertex <= size; vertex++ ) {
        _vertices[vertex-1] = vertex;
        _edges[vertex].reserve(dimacs_graph.degrees[vertex]);
    }

    for ( const std::pair<int, int>& edge : dimacs_graph.edges ) {
        // skipping already inserted edges
        if ( std::find(_edges[edge.first].begin(), 
                       _edges[edge.first].end(), 
                        edge.second) != _edges[edge.first].end() ) {
            continue;
        }
        if ( edge.first != edge.second ) {
            _edges[edge.first].push_back(edge.second);
            _edges[edge.second].push_back(edge.first);
        } else if ( edge.first == edge.second ) {
            _edges[edge.first].push_back(edge.second);
        } 
    }

    _degrees.clear();
    _degrees.resize(_vertices.size()+1);
    for (int i = 0; i < _vertices.size(); i++) {
        _degrees[_vertices[i]] = _edges[_vertices[i]].size();
    }

}

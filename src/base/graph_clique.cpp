#include "graph_clique.hpp"

GraphClique::GraphClique(int num_vertices) : adjacency_list_1(num_vertices) {}

void GraphClique::AddEdge_B(int u, int v) {
    adjacency_list_1[u].insert(v);
    adjacency_list_1[v].insert(u);
}

int GraphClique::GetNumVertices() const {
    return adjacency_list_1.size();
}

const std::set<int>& GraphClique::GetNeighbors_B(int v) const {
    return adjacency_list_1[v];
}

// Constructor: Initializes an empty graph
GraphClique::GraphClique() {}

// Adds an edge between two vertices with an optional weight (default is 1)
void GraphClique::AddEdge(int u, int v, int weight) {
    adjacency_list_[u].insert(v);
    adjacency_list_[v].insert(u);

    // Assign default weight if not already set
    if (vertex_weights_.find(u) == vertex_weights_.end()) {
        vertex_weights_[u] = weight;
    }
    if (vertex_weights_.find(v) == vertex_weights_.end()) {
        vertex_weights_[v] = weight;
    }
}

// Returns the set of all vertices in the graph
std::set<int> GraphClique::GetVertices() const {
    std::set<int> vertices;
    for (const auto& entry : adjacency_list_) {
        vertices.insert(entry.first);
    }
    return vertices;
}

// Returns the neighbors of a given vertex
std::set<int> GraphClique::GetNeighbors(int v) const {
    if (adjacency_list_.find(v) != adjacency_list_.end()) {
        return adjacency_list_.at(v);
    }
    return {};
}

// Checks if two vertices are adjacent
bool GraphClique::IsNeighbor(int u, int v) const {
    if (adjacency_list_.find(u) != adjacency_list_.end()) {
        return adjacency_list_.at(u).count(v) > 0;
    }
    return false;
}

// Returns the weight of a given vertex (default is the sum of its edges)
int GraphClique::GetVertexWeight(int v) const {
    if (vertex_weights_.find(v) != vertex_weights_.end()) {
        return vertex_weights_.at(v);
    }
    return 0;
}

// Returns the sum of weights of all neighbors of a vertex
int GraphClique::GetNeighborsWeightSum(int v) const {
    int sum = 0;
    if (adjacency_list_.find(v) != adjacency_list_.end()) {
        for (int neighbor : adjacency_list_.at(v)) {
            sum += GetVertexWeight(neighbor);
        }
    }
    return sum;
}

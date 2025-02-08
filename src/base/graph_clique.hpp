#ifndef GRAPH_CLIQUE_HPP
#define GRAPH_CLIQUE_HPP

#include <vector>
#include <set>
#include <unordered_map>

class GraphClique {
public:
    explicit GraphClique(int num_vertices);
    void AddEdge_B(int u, int v);
    int GetNumVertices() const;
    const std::set<int>& GetNeighbors_B(int v) const;
    // Constructor: Initializes an empty graph
    GraphClique();

    // Adds an edge between two vertices with optional weight (default: 1)
    void AddEdge(int u, int v, int weight = 1);

    // Returns the set of all vertices in the graph
    std::set<int> GetVertices() const;

    // Returns the neighbors of a given vertex
    std::set<int> GetNeighbors(int v) const;

    // Returns whether two vertices are adjacent
    bool IsNeighbor(int u, int v) const;

    // Returns the weight of a given vertex (default: sum of its edges)
    int GetVertexWeight(int v) const;

    // Returns the sum of weights of all neighbors of a vertex
    int GetNeighborsWeightSum(int v) const;


private:
    std::vector<std::set<int> > adjacency_list_1;
    std::unordered_map<int, std::set<int>> adjacency_list_;  // Stores adjacency list
    std::unordered_map<int, int> vertex_weights_;  // Stores vertex weights

};

#endif // GRAPH_CLIQUE_HPP

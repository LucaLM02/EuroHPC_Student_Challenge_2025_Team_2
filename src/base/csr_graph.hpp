#ifndef CSR_GRAPH_HPP
#define CSR_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

class CSRGraph : public Graph {

private:
    int nEdges;
    int nVertices;
    std::vector<int> edges;
    std::vector<int> offsets;

public:
    CSRGraph(const Dimacs& dimacs_graph);
    CSRGraph(const CSRGraph& other);

    void AddEdge(int v, int w) override;
    void RemoveEdge(int v, int w) override;
    void AddVertex(int v) override;
    void RemoveVertex(int v) override;
    bool MergeVertices(int v, int w) override;

    void GetNeighbours(int vertex, std::vector<int> &result) const override;
    void GetNeighbours(int vertex, std::set<int> &result) const override;
    int GetNeighboursIndex(int vertex) const;

    bool HasEdge(int v, int w) const = 0;
    void GetEdges(int vertex, std::set<std::pair<int,int>> &result) const override;
    void GetEdges(int vertex, std::vector<std::pair<int, int>> &result) const override;

    size_t GetNumVertices() const override;
    size_t GetNumEdges() const override;

};

#endif // CSR_GRAPH_HPP
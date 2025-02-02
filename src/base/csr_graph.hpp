#ifndef CSR_GRAPH_HPP
#define CSR_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

/*
    TODO: handle properly vertex removal: since offset is accessed as offset[vertex], if
          vertices aren't renamed when vertex is removed, then the entry of `vertex` must
          remain
    TODO: otherwise could be possible to impose a rewriting of the vertices when the topology changes
        - like RemoveVertex(int w, bool reorder=false);
    TODO: write docs
    TODO: vertices numbering convention starts from 1, so an extra space at index 0 is added but not used
*/

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
    void MergeVertices(int v, int w) override;

    void GetNeighbours(int vertex, std::vector<int> &result) const override;
    void GetNeighbours(int vertex, std::set<int> &result) const override;
    int GetNeighboursIndex(int vertex) const;

    bool HasEdge(int v, int w) const = 0;

    size_t GetNumVertices() const override;
    size_t GetNumEdges() const override;

};

#endif // CSR_GRAPH_HPP
#ifndef CSR_GRAPH_HPP
#define CSR_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

class CSRGraph : public Graph {

private:
    std::vector<int> edges;
    std::vector<int> offsets;

public:
    CSRGraph(const Dimacs& dimacs_graph);
    CSRGraph(const CSRGraph& other);
    void addEdge(int v, int w) override;
    std::vector<int> getNeighbours(int vertex) const override;
    size_t getNumVertices() const override;
    size_t getNumEdges() const override;
    bool merge(int v, int w) override;

    int getNeighboursIndex(int vertex) const;
};

#endif // CSR_GRAPH_HPP
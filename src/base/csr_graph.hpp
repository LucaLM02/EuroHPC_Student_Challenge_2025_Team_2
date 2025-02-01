#ifndef CSR_GRAPH_HPP
#define CSR_GRAPH_HPP

#include "graph.hpp"

class CSRGraph : public Graph {

private:
    std::vector<int> edges;
    std::vector<int> offsets;

public:
    CSRGraph();
    void addEdge(int v, int w) override;
    std::vector<int> getNeighbours(int vertex) const override;
    size_t getNumVertices() const override;
    size_t getNumEdges() const override;
    bool merge(int v, int w) override;
    ~Graph() =default;

    int getNeighboursIndex(int vertex) const;
};

#endif // CSR_GRAPH_HPP
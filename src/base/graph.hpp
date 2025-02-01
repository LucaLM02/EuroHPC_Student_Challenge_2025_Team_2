#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>

Class Graph {
public:
    Graph();
    virtual void addEdge(int v, int w) =0;
    virtual std::vector<size_t> getNeighbours(int vertex) const = 0;
    virtual size_t getNumVertices() const = 0;
    virtual size_t getNumEdges() const = 0;
    virtual bool merge(int v, int w) = 0;
    virtual ~Graph();
private:
    int nEdges;
    int nVertices;
};

#endif // GRAPH_HPP
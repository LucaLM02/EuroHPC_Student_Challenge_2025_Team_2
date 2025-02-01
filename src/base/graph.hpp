#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>

class Graph {
    public:
        Graph() {};
        virtual void addEdge(int v, int w) = 0;
        virtual std::vector<int> getNeighbours(int vertex) const = 0;
        virtual size_t getNumVertices() const = 0;
        virtual size_t getNumEdges() const = 0;
        virtual bool merge(int v, int w) = 0;
        virtual ~Graph() = default;

    protected:
        int nEdges;
        int nVertices;
};

#endif // GRAPH_HPP
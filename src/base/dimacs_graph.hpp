#ifndef DIMACS_GRAPH_HPP
#define DIMACS_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

/*
    TODO: handle the fact that only symmetric vertices are
*/

class DimacsGraph : public Graph {
    public:
        DimacsGraph(Dimacs& dimacs);

        // -------------------- MODIFIERS --------------------
        virtual void AddEdge(int v, int w) override;
        virtual void RemoveEdge(int v, int w) override;
        virtual void AddVertex(int v) override;
        virtual void RemoveVertex(int v) override;
        virtual void MergeVertices(int v, int w) override;

        // --------------------- GETTERS ----------------------
        virtual void GetNeighbours(int vertex, std::vector<int> &result) const override;
        virtual void GetNeighbours(int vertex, std::set<int> &result) const override;

        virtual bool HasEdge(int v, int w) const override;
        virtual void GetVertices(std::set<int> &result) const;
        virtual void GetVertices(std::vector<int> &result) const;

        virtual size_t GetNumVertices() const override;
        virtual size_t GetNumEdges() const override;
        virtual ~DimacsGraph() = default;

    public:
        Dimacs&       _dimacs;
        std::set<int> _vertices;

};

#endif // DIMACS_GRAPH_HPP

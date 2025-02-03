#ifndef DIMACS_GRAPH_HPP
#define DIMACS_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

class DimacsGraph : public Graph {
    public:
        DimacsGraph(Dimacs& dimacs);

        // -------------------- MODIFIERS --------------------
        virtual void AddEdge(int v, int w) override;
        virtual void RemoveEdge(int v, int w) override;
        virtual void AddVertex(int v) override;
        virtual void RemoveVertex(int v) override;
        virtual void RemoveVertexWithRenaming(int v) override;
        virtual void MergeVertices(int v, int w) override;

        // --------------------- GETTERS ----------------------
        virtual void GetNeighbours(int vertex, std::vector<int> &result) const override;
        virtual void GetNeighbours(int vertex, std::set<int> &result) const override;

        virtual bool HasEdge(int v, int w) const override;
        virtual void GetUnorderedVertices(std::set<int> &result) const override;
        virtual const std::vector<int>& GetVertices() const override;

        virtual void SetVertices(std::vector<int>& vertices);

        virtual size_t GetNumVertices() const override;
        virtual size_t GetNumEdges() const override;

        virtual unsigned int GetDegree(int vertex) const;
        virtual std::vector<int> GetDegrees() const;
        virtual unsigned int GetMaxDegree() const;
        virtual int GetVertexWithMaxDegree() const;

        virtual ~DimacsGraph() = default;

    private:
        /**
         * @brief removes the vertex without touching offset, degrees and edges
         * 
         * @param vertex 
         */
        void _RemoveOnlyVertex(int vertex);
        void _RemoveDegree(int vertex);

        Dimacs&       _dimacs;
        std::vector<int> _vertices;

};

#endif // DIMACS_GRAPH_HPP

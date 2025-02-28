#ifndef DIMACS_GRAPH_HPP
#define DIMACS_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

#include <memory>

class DimacsGraph : public Graph {
    public:

        static DimacsGraph* LoadFromDimacs(const std::string& file_name);
        DimacsGraph() = default;

        bool isEqual(const Graph &ot) const override{
            return false;
        }
        // -------------------- MODIFIERS --------------------
        
        virtual void AddHistory(GraphHistory graph_history) {};
        virtual void AddEdge(int v, int w) override;
        virtual void RemoveEdge(int v, int w) override;

        virtual int AddVertex() override;
        virtual void RemoveVertex(int v) override;
        virtual void SetVertices(std::vector<int>& vertices);

        virtual void MergeVertices(int v, int w) override;

        virtual void SetColoring(const std::vector<unsigned short>& colors) override {};
        virtual void SetColoring(int vertex, unsigned short color) override {};
        virtual void SetFullColoring(const std::vector<unsigned short>& colors) override {}
        virtual void ClearColoring() override {};

        // -------------------- ORDERING ----------------------

        virtual void SortByDegree(bool ascending=true) override {}
        virtual void SortByExDegree(bool ascending=true) override {}
        virtual void SortByColor(bool ascending=true) override {};


        // --------------------- GETTERS ----------------------
        virtual void GetNeighbours(int vertex, std::vector<int> &result) const override;
        virtual void GetNeighbours(int vertex, std::set<int> &result) const override;

        virtual bool HasEdge(int v, int w) const override;

        virtual void GetUnorderedVertices(std::set<int> &result) const override;
        virtual const std::vector<int>& GetVertices() const override;
        virtual int GetVertexByIndex(int index) const override;
        virtual int GetHighestVertex() const override;

        virtual size_t GetNumVertices() const override;
        virtual size_t GetNumEdges() const override;

        virtual unsigned int GetDegree(int vertex) const override;
        virtual std::vector<int> GetDegrees() const override;
        virtual void GetDegrees(std::vector<int>& result) const override;
        virtual std::vector<int> GetFullDegrees() const override { return {}; };
        virtual void GetFullDegrees(std::vector<int>& result) const override {};
        virtual unsigned int GetMaxDegree() const override;
        virtual int GetVertexWithMaxDegree() const override;
        virtual int GetExDegree(int vertex) const override { return 0;};

        virtual std::vector<int> GetMergedVertices(int vertex) const override { return {}; };

        virtual std::vector<unsigned short> GetColoring() const override { return {}; };
        virtual std::vector<unsigned short> GetFullColoring() const override { return {}; };
        virtual unsigned short GetColor(int vertex) const override { return 0; };

        // -------------------- SERIALIZATION --------------------
        std::string Serialize() const;
        void Deserialize(const std::string& data);

        virtual std::unique_ptr<Graph> Clone() const override;

        virtual ~DimacsGraph() = default;

    protected:
        DimacsGraph(const std::string& file_name);

    private:
        /**
         * @brief removes the vertex without touching offset, degrees and edges
         * 
         * @param vertex 
         */
        void _RemoveOnlyVertex(int vertex);
        void _RemoveDegree(int vertex);

        Dimacs       _dimacs;
        std::vector<int> _vertices;
        std::set<int> _deleted_vertices;
        mutable std::vector<int> _tmp_degrees;

};

#endif // DIMACS_GRAPH_HPP

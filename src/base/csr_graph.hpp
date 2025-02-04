#ifndef CSR_GRAPH_HPP
#define CSR_GRAPH_HPP

#include "graph.hpp"
#include "dimacs.hpp"

/*
    TODO: otherwise could be possible to impose a rewriting of the vertices when the topology changes
        - like RemoveVertex(int w, bool reorder=false);
    TODO: write docs
*/

class CSRGraph : public Graph {
    public:
        static CSRGraph* LoadFromDimacs(const std::string& file_name);

        CSRGraph()=default;
        CSRGraph(const CSRGraph& other)=default;

        // -------------------- MODIFIERS --------------------
        virtual void AddEdge(int v, int w) override;
        virtual void RemoveEdge(int v, int w) override;

        virtual int AddVertex() override;
        virtual void RemoveVertex(int v) override;
        virtual void RemoveVertexWithRenaming(int v) override;
        virtual void SetVertices(std::vector<int>& vertices);

        virtual void MergeVertices(int v, int w) override;

        /**
         * @brief orders the vertices by degree. 
         * 
         * @param ascending 
         */
        virtual void OrderByDegree(bool ascending=true) {};
        /**
         * @brief recompress the graph by removing data which has been deleted but not 
         * effectively erased not to decrease performance
         * @note this method should be called when ~100 vertices are removed from the graph. 
         * Too frequent invocations of this method aren't efficient
         * @details when removing a node (due to RemoveVertex or OrderByDegree), some spaces 
         * are left in the internal memory representation. This avoids burdersome calculations 
         * but if the vertices removed are a lot, the spaces accumulate and make this graph 
         * representation less efficient.
         */
        virtual void Compress() {};

        // --------------------- GETTERS ----------------------
        virtual void GetNeighbours(int vertex, std::vector<int> &result) const override;
        virtual void GetNeighbours(int vertex, std::set<int> &result) const override;

        virtual bool HasEdge(int v, int w) const override;

        virtual void GetUnorderedVertices(std::set<int> &result) const override;
        virtual const std::vector<int>& GetVertices() const override;
        virtual int GetVertexByIndex(int index) const;
        virtual const std::set<int>& GetDeletedVertices() const override;

        virtual size_t GetNumVertices() const override;
        virtual size_t GetNumEdges() const override;

        virtual unsigned int GetDegree(int vertex) const;
        virtual const std::vector<int>& GetDegrees() const;
        virtual void GetDegrees(std::vector<int>& result) const;
        virtual unsigned int GetMaxDegree() const;
        virtual int GetVertexWithMaxDegree() const;

        virtual ~CSRGraph() = default;

    private:
        CSRGraph(const Dimacs& dimacs_graph);
        void _ComputeCacheDegrees() const;
        size_t _nEdges;
        std::vector<int> _vertices;         // TODO: fill this
        std::set<int> _removed_vertices;
        std::vector<std::vector<int>> _edges;
        mutable std::vector<int> _cache_degrees;
        mutable bool _cache_degrees_invalidated = true;

        int              _index_of_max_degree;


};

#endif // CSR_GRAPH_HPP
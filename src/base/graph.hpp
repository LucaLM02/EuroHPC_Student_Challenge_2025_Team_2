#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <set>

/**
 *  @brief Abstract class that represents an undirected (possibly loop-)graph, composed by a set of vertices and edges. 
 * 
 *  @details
 *  Simple abstract class that represents a graph, composed by vertices, which can assume 
 *  any value > 0, and edges, represented (from user perspective) as a std::pair<int, int>
 *  <br>
 *  The graph can be modified, since it is useful for applying Zykov coloring algorithm.
 *  Vertex and edges can be added and removed. 2 vertices can be merged (see MergeVertices 
 *  for a more detailed description)
 */
class Graph {
    public:
        /**
         * @brief Destroys the Graph object
         */
        virtual ~Graph() = default;

        /**
         *  @brief adds the edge <v,w>=<w,v> to the graph
         * 
         *  @warning for efficiency reasons, it is user responsability to ensure
         *           that the edge wasn't already present. If this is not true,
         *           the edge will be duplicated, leading possibly to undesired
         *           behaviours
         */
        virtual void AddEdge(int v, int w) = 0;

        /**
         *  @brief removes the edge <v,w>=<w,v> from the graph
        */
        /**
         *  @brief removes the edge <v,w>=<w,v> from the graph
         */
        virtual void RemoveEdge(int v, int w) = 0;
        /**
         *  @brief adds a vertex to the graph
         * 
         *  @warning for efficiency reasons, it is user responsability to ensure
         *           that the vertex wasn't already present. If this is not true,
         *           the vertex will be duplicated, leading possibly to undesired
         *           behaviours
         */
        virtual void AddVertex(int v) = 0;
        /**
         *  @brief removes a vertex from the graph without changing the other vertices.
         * 
         *  @details
         *  Removes a vertex from the graph without changing the other vertices. <br>
         *  In particular, it DOES NOT change the names/tags of other vertices
         */
        virtual void RemoveVertex(int v) = 0;

        /**
         *  @brief removes a vertex from the graph renaming other vertices
         * 
         *  @details
         *  Removes a vertex from the graph. <br>
         *  Opposed to RemoveVertex, it renames the other vertices in such a way that the final
         *  vertices are contiguous from 1 to GetNumVertices() (both included).
         */
        virtual void RemoveVertexWithRenaming(int v) = 0;

        /**
         *  @brief merges vertex v and vertex w into one single vertex, which will be called vertex v
         *  @details
         *  Merges vertex v and vertex w into one single vertex, which will be called vertex v
         *  Merging means that the final vertex will have all the neighbours of v and all 
         *  the neighbours of w. <br>
         *  @note If the 2 vertices were neighbours, the final vertex `v` will have a loop
         */
        virtual void MergeVertices(int v, int w) = 0;

        /**
         *  @brief gets the neighbours of `vertex` as a vector
         *  @details
         *  Gets the neighbours of `vertex` as a vector <br>
         *  Child classes might impose a certain order in the result
         */
        virtual void GetNeighbours(int vertex, std::vector<int> &result) const = 0;
        /**
         *  @brief gets the neighbours of `vertex` as a set
         */
        virtual void GetNeighbours(int vertex, std::set<int> &result) const = 0;

        /**
         * @brief returns true iff <v,w>=<w,v> is an edge of this graph
         * 
         * @param v first vertex
         * @param w second vertex
         * @return true if <v,w>=<w,v> edge belongs to this graph
         * @return false if <v,w>=<w,v> edge doesn't belong to this graph
         */
        virtual bool HasEdge(int v, int w) const = 0;

        /**
         *  @brief returns set of vertices used in the graph but without a precise order
         *  @details
         *  Returns set of vertices used in the graph without a precise order<br>
         *  Note that vertices aren't necessarily the number between 0 and GetNumVertices(),
         *  they shall be any sequence of numbers >= 0
         */
        virtual void GetUnorderedVertices(std::set<int> &result) const = 0;
        /**
         *  @brief returns vector of vertices used in the graph. The order matters here
         *  @details
         *  Returns vector of vertices used in the graph <br>
         *  Note that vertices aren't necessarily the number between 0 and GetNumVertices(),
         *  they shall be any sequence of numbers >= 0
         */
        virtual const std::vector<int>& GetVertices() const = 0;

        /**
         *  @brief sets the vertices of the graph. Used typically to set the order of them
         *  @warning for efficiency reasons, this method does not delete the edges which do not have one of the 2 vertices after this method is called
         */
        virtual void SetVertices(std::vector<int>& vertices) = 0;

        virtual size_t GetNumVertices() const = 0;
        virtual size_t GetNumEdges() const = 0;
        /**
         * @brief Get the degree of `vertex`
         * 
         * @param vertex vertex of which the degree is returned
         * @return unsigned int degree of `vertex`
         */
        virtual unsigned int GetDegree(int vertex) const = 0;
        /**
         * @brief Returns an array of degrees, one for each vertex
         * Note: no order is imposed
         * @return std::vector<int> array of degrees
         */
        virtual std::vector<int> GetDegrees() const = 0;
        /**
         * @brief Returns the maximum degree (delta(G)) among the degree of all vertices of this graph
         * 
         * @return unsigned int maximum degree of the graph (delta(G))
         */
        virtual unsigned int GetMaxDegree() const = 0;
        /**
         * @brief Gets the vertex with the max degree 
         * 
         * @return int the vertex with the max degree
         */
        virtual int GetVertexWithMaxDegree() const = 0;

};

#endif // GRAPH_HPP

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
         *  @brief returns set of vertices used in the graphc

         *  @details
         *  Returns set of vertices used in the graph <br>
         *  Note that vertices aren't necessarily the number between 0 and GetNumVertices(),
         *  they shall be any sequence of numbers >= 0
         */
        virtual void GetVertices(std::set<int> &result) const = 0;
        /**
         *  @brief returns vector of vertices used in the graph

         *  @details
         *  Returns vector of vertices used in the graph <br>
         *  Note that vertices aren't necessarily the number between 0 and GetNumVertices(),
         *  they shall be any sequence of numbers >= 0
         */
        virtual void GetVertices(std::vector<int> &result) const = 0;

        virtual size_t GetNumVertices() const = 0;
        virtual size_t GetNumEdges() const = 0;


};

#endif // GRAPH_HPP
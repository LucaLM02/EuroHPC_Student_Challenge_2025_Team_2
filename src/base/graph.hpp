#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <set>
#include <memory>

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
 * 
 *  @warning for efficiency reasons, children classes might not be checking the correctness
 *  of the input parameter. See each method for a more specific warning
 */
class Graph {
    public:
        /**
         * @brief Destroys the Graph object
         */
        virtual ~Graph() = default;

        virtual bool isEqual(const Graph &ot) const = 0;

        // ================================== MODIFIERS ==================================
        /**
         *  @brief adds the edge <v,w>=<w,v> to the graph
         * 
         *  @warning for efficiency reasons, it is user responsability to ensure
         *           that the edge wasn't already present and that v,w>0. 
         *           If this is not true, the class behaviour is undefined
         */
        virtual void AddEdge(int v, int w) = 0;

        /**
         *  @brief removes the edge <v,w>=<w,v> from the graph
         *  @warning undefined behaviour if either v or w don't belong to this graph vertices
         */
        virtual void RemoveEdge(int v, int w) = 0;


        /**
         * @brief sets the vertices of the graph. Used typically to set the order of them
         * @warning undefined behaviour if vertices is not a permutation of 
         *          this->GetVertices()
         */
        virtual void SetVertices(std::vector<int>& vertices) = 0;

        /**
         *  @brief adds a vertex to the graph
         * 
         */
        virtual int AddVertex() = 0;
        
        /**
         *  @brief removes a vertex from the graph without changing the other vertices.
         * 
         *  @details
         *  Removes a vertex from the graph without changing the other vertices. <br>
         *  In particular, it DOES NOT change the names/tags of other vertices

         *  @warning for efficiency reasons, does not check whether v and w are
         *           part of the vertices of the graph. If v (and or w) doesn't
         *           belong to the vertices of this graph, the behavious is not
         *           determined
         */
        virtual void RemoveVertex(int v) = 0;

        /**
         *  @brief merges vertex v and vertex w into one single vertex, which will be called vertex v
         *  @details
         *  Merges vertex v and vertex w into one single vertex, which will be called vertex v
         *  Merging means that the final vertex will have all the neighbours of v and all 
         *  the neighbours of w. <br>
         *  @note If the 2 vertices were neighbours, the final vertex `v` will have a loop
         *  @warning undefined behaviour if either v or w do not belong to this graph vertices
         */
        virtual void MergeVertices(int v, int w) = 0;


        /**
         * @brief sets the coloring of the graph, i.e. one color for each vertex
         * 
         * @param colors color vector, which should be ordered as vertices are. 
         *               Each color must be >= 0. Color = 0 means that no color is
         *               assigned to that vertex
         */
        virtual void SetColoring(const std::vector<unsigned short>& colors) = 0;

        /**
         * @brief sets the color of a single vertex
         * 
         * @param vertex vertex to which color is being set
         * @param color color to set (>= 0, if = 0 then it has no color)
         */
        virtual void SetColoring(int vertex, unsigned short color) = 0;

        /**
         * @brief sets the coloring of the graph, i.e. one color for each vertex. 
         *        The difference arieses in how colors vector is composed.
         *        If the graph has non continuous numbering, `colors` will have
         *        entries set to 0 which aren't used.
         *        However, this method is typically computationally ligher for
         *        both the user of this class and this class
         * 
         * @param colors color vector. For each vertex its color is accessible as
         *               colors[vertex]. This means that if the graph has a non 
         *               contiguous numbering, `colors` will have entries set to 0.
         */
        virtual void SetFullColoring(const std::vector<unsigned short>& colors) = 0;

        /**
         * @brief clears the coloring of the graph. 
         */
        virtual void ClearColoring() = 0;

        // ================================== ORDERING ===================================

        //virtual void Order(GraphOrdering ...)

        /**
         * @brief orders the vertices by degree. 
         * 
         * @param ascending ascending or descending order
         */
        virtual void SortByDegree(bool ascending=false) = 0;

        /**
         * @brief orders the vertices by the ex degree (see long description)
         * @details ex-degree = vertex degree + sum of the degrees of the neighbours
         * @param ascending ascending or descending order
         */
        virtual void SortByExDegree(bool ascending=false) = 0;

        /**
         * @brief orders the vertices by their color
         * 
         * @param ascending ascending or descending order
         */
        virtual void SortByColor(bool ascending=false) = 0;


        // ================================== GETTERS ====================================
        /**
         *  @brief gets the neighbours of `vertex` as a vector
         *  @details
         *  Gets the neighbours of `vertex` as a vector <br>
         *  Child classes might impose a certain order in the result
         *  @warning undefined behaviour if v doesn't belong to this graph vertices
         */
        virtual void GetNeighbours(int vertex, std::vector<int> &result) const = 0;
        /**
         *  @brief gets the neighbours of `vertex` as a set
         *  @warning undefined behaviour if v doesn't belong to this graph vertices
         */
        virtual void GetNeighbours(int vertex, std::set<int> &result) const = 0;

        /**
         * @brief returns true iff <v,w>=<w,v> is an edge of this graph
         * 
         * @param v first vertex
         * @param w second vertex
         * @return true if <v,w>=<w,v> edge belongs to this graph
         * @return false if <v,w>=<w,v> edge doesn't belong to this graph
         * @warning undefined behaviour if either v or w don't belong to this graph vertices
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
         *  @brief returns vector of vertices used in the graph. The order matters here.
         *         When a vertex is added or deleted, the result might be updated 
         *         "automatically" since it is a reference to an inner object
         *  @details
         *  Returns vector of vertices used in the graph <br>
         *  Note that vertices aren't necessarily the number between 0 and GetNumVertices(),
         *  they shall be any sequence of numbers >= 0
         */
        virtual const std::vector<int>& GetVertices() const = 0;

        /**
         * @brief with vectors such as the one returned by GetDegrees(), user might want to
         *        know each element (so at each index) which vertex is assigned.
         *        This is the method for gathering that information
         * @param index index to map to vertex
         * @return the vertex associated to index
         * @warning the behaviour is undefined if index is not associated to any vertex
         */
        virtual int GetVertexByIndex(int index) const = 0;

        /**
         * @brief returns the vertex with the highest value in the graph
         * 
         * @return int 
         */
        virtual int GetHighestVertex() const = 0;

        /**
         * @brief returns the vertices deleted from the graph. When a vertex is deleted,
         *        the result is "automatically" updated, since it is a reference to an
         *        inner object
         * @note when the graph vertices are renamed with values from 1 to GetNumVertices(), 
         *       then this this function will return an empty set, until a new vertex is 
         *       deleted
         * @return std::set<int>
         */
        virtual const std::set<int>& GetDeletedVertices() const = 0;

        virtual size_t GetNumVertices() const = 0;
        virtual size_t GetNumEdges() const = 0;
        /**
         * @brief Get the degree of `vertex`
         * 
         * @param vertex vertex of which the degree is returned
         * @return unsigned int degree of `vertex`
         * @warning undefined behaviour if v doesn't belong to this graph vertices
         */
        virtual unsigned int GetDegree(int vertex) const = 0;
        /**
         * @brief Returns an array of degrees, ordered as the vertices are
         * @note Also, don't to degrees[vertex], see GetFullDegrees() for this
         * @return std::vector<int> array of degrees
         * @see GetVertexByIndex for understanding at which vertex a certain degree 
         *      is associated
         */
        virtual std::vector<int> GetDegrees() const = 0;
        /**
         * @brief returns, through the reference parameter, the degrees in the same order 
         *        of the vertices
         * @param result the degree's vector to fill
         * @see GetVertexByIndex for understanding at which vertex a certain degree 
         *      is associated
         */
        virtual void GetDegrees(std::vector<int>& result) const = 0;
        /**
         * @brief returns, through the reference parameter, the degrees such 
         *        that result[vertex] -> vertex_degree
         * @return std::vector<int> array of degrees
         */
        virtual std::vector<int> GetFullDegrees() const = 0;
        /**
         * @brief returns, through the reference parameter, the degrees such 
         *        that result[vertex] -> vertex_degree
         * @param result the degree's vector to fill
         */
        virtual void GetFullDegrees(std::vector<int>& result) const = 0;

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
        /**
         * @brief gets the sum of the degrees of the vertex neighbourhood
         * 
         * @return int the sum of the degrees of the vertex neighbourhood
         */
        virtual int GetExDegree(int vertex) const = 0;

        /**
         * @brief Get the vertices that were merged into the given vertex. The first element is the vertex itself
         *  
         * @param vertex 
         */
        virtual std::vector<int> GetMergedVertices(int vertex) const = 0;

        /**
         * @brief gets the coloring of the graph
         * 
         * @return std::vector<unsigned short> vector of colors, ordered as the vertices are
         */
        virtual std::vector<unsigned short> GetColoring() const = 0;

        /**
         * @brief gets the coloring of the graph, i.e. one color for each vertex. 
         *        The difference arieses in how colors vector is composed.
         *        If the graph has non continuous numbering, `colors` will have
         *        entries set to 0 which aren't used.
         *        However, this method is typically computationally ligher for
         *        both the user of this class and this class
         * 
         * @returns The color vector. For each vertex its color is accessible as
         *          colors[vertex]. This means that if the graph has a non 
         *          contiguous numbering, `colors` will have entries set to 0.
         */
        virtual std::vector<unsigned short> GetFullColoring() const = 0;
        
        /**
         * @brief gets the color of a certain vertex. If vertex is not in the graph or
         *        if the vertex is not colored, 0 is returned, otherwise a value > 0.
         * 
         * @param vertex vertex of which the color is returned
         * @return unsigned short color of the vertex
         */
        virtual unsigned short GetColor(int vertex) const = 0;

        // ================================= COPYING ====================================
        virtual std::unique_ptr<Graph> Clone() const = 0;
        
        // ================================= SERIALIZATION ==============================
        virtual std::string Serialize() const = 0;
        virtual void Deserialize(const std::string& data) = 0;
};

#endif // GRAPH_HPP

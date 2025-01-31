#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H


#include <vector>
#include <string>


/**
 * This file defines a generic graph loader class (interface)
 * */
 
 
 /**
  * @class GraphLoader
  * @author Matjaž
  * @date 11/11/16
  * @file GraphLoader.h
  * @brief The abstract class GraphLoader which all loaders should extend.
  */
 class GraphLoader {
public:
    virtual ~GraphLoader() {};
    /**
     * @brief Load a file specified by @link "filename.
     * On errors, false is returned and the details of error are stored and can be queried by function getError();
     * @param fname the file name (and path)
     * @return true if loads ok, false if not (in that case see @link #getError() for details)
     */
    virtual bool load(const char* fname) = 0;
    
    /**
     * @brief get the number of vertices in the graph; graph must be loaded first or the result will be undefined.
     * If the graph was not loaded, the return will be undefined.
     * @return the total number of vertices
     */
    virtual unsigned int getNumVertices() const = 0;
    
    /**
     * @brief get the maximum vertex index; usually this will be number of vertices - 1, but it could also be number of vertices for 1-based indexing.
     * If the graph was not loaded, the return will be undefined.
     * @return the maximum index that is represents a valid vertex
     */
    virtual unsigned int getMaxVertexIndex() const = 0;
    
    /**
     * @brief the number of edges in a graph (undirected graphs should count each edge twice: from _a_ to _b_ and from _b_ to _a_)
     * If the graph was not loaded the will be undefined.
     * @return the number of edges (directional)
     */
    virtual unsigned int getNumEdges() const = 0;
    
    /**
     * @brief get (calculate if necessary) the adjacency matrix, where each edge is represented by and 8bit boolean variable.
     * If the graph was not loaded this should return an empty (0x0) matrix.
     * @return the adjacency matrix
     */
    virtual std::vector<std::vector<char> > getAdjacencyMatrix() const = 0;
    
    /**
     * @brief Get the vertex degrees vector.
     * If the graphs was not loaded, the vector will be empty.
     * @return the vector of vertex degrees.
     */
    virtual std::vector<int> getDegrees() const = 0;
    
    /**
     * @brief Calculate some graph statistics: min and max degree, and histogram of degrees.
     * @param maxdegree the output variable - minimal degree of a vertex in the graph
     * @param minDegree the output variable - maximal degree of a vertex in the graph
     * @param degreeHistogram - histogram of degrees, one entry for each degree between min and max (inclusive).
     */
    virtual void calculateGraphStats(int& maxdegree, int& minDegree, std::vector<float>& degreeHistogram);
    
    /**
     * @brief Get the last error message (load produces error messages)
     * @return the last error message
     */
	virtual const std::string& getError() const = 0;
    
    /**
     * @brief Will return true if files contain 1-based indexing of vertices and this is converted to 0-based indexing on load
     * @return true for 1-based indexing converted to 0-based indexing
     */
    virtual bool verticesAreMappedFrom1based() const = 0;
};


#endif //GRAPH_LOADER_H

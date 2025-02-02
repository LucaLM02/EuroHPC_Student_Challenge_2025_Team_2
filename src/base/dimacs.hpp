#ifndef DIMACS_H
#define DIMACS_H


#include <utility>
#include <vector>
#include <string>
#include "graph_loader.hpp"


class Dimacs : public GraphLoader {
public:
    typedef std::pair<int, int> Edge;
    
protected: public:
    //std::vector<Edge> edges;
    std::vector<Edge> edges;
    
    std::vector<int> degrees;
    unsigned int numVertices;
    unsigned long maxVertexIndex;
    unsigned long long int adjacencyMatrixSizeLimit; // TODO: set this from commandline
	std::string error;
	bool errorFlag;
	bool edgeNotSpecified = false;
    
public:
    Dimacs();
    ~Dimacs();
    
    void allowNotSpecifiedEdge(bool allow=true) {edgeNotSpecified = allow;}
    bool load(const char* fname) override;
    unsigned int getNumVertices() const override        {return numVertices;}
    unsigned int getMaxVertexIndex() const override    {return maxVertexIndex;}
    unsigned int getNumEdges() const override           {return edges.size();}
    std::vector<std::vector<char> > getAdjacencyMatrix() const override;
    std::vector<int> getDegrees() const override        {return degrees;}
	const std::string& getError() const override        {return error;}
    bool verticesAreMappedFrom1based() const override   {return false;}
	
private:
	bool parseProblemLine(std::istringstream& ss, std::string& line);
	bool parseSpecsLine(std::istringstream& ss, std::string& line);
    bool parseLine(std::istringstream& ss, std::string& line);
};


#endif // DIMACS_H


#include "dimacs.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>


const bool ENABLE_DEBUG_OUT_READING = false;


Dimacs::Dimacs() : numVertices(0), adjacencyMatrixSizeLimit(1000000000 /* 1GB */) {
}

Dimacs::~Dimacs() {
}

bool Dimacs::load(const char* fname) {
    std::ifstream file(fname);

    if (!file) {
        error = "ifstream invalid - file cannot be opened for reading";
        return false;
    }
	
    //std::cerr << "loading DIMACS, ..." << std::endl;
	errorFlag = false;
	
	// parse function should return true when it is done parsing, false if it is not done yet
	typedef bool (Dimacs::*parseFunc)(std::istringstream&, std::string&);
	// loop function will return true on success, false otherwise 
	auto loop = [&file, this](parseFunc f)->bool {
		std::string line;
		while (file) {
			if (!std::getline(file, line))
                break;
            if (ENABLE_DEBUG_OUT_READING)
                std::cout << "getline " << line << "\n";
			if (line.size() > 0) {
				std::istringstream ss(line);
				ss.ignore(2);
				if (!(this->*f)(ss, line)) {
                    error = "parse function returned error on line '" + line + "': \"" + error + "\"";
					return false;                    
                }
			}
		}
		return !errorFlag;
	};
	
	return loop(&Dimacs::parseLine);
}


std::vector<std::vector<char> > Dimacs::getAdjacencyMatrix() const {
    if (getMaxVertexIndex() * getMaxVertexIndex() > adjacencyMatrixSizeLimit)
        throw std::runtime_error("Cannot create adjacency matrix because the number of vertices is to large");
    std::vector<std::vector<char> > matrix;
    matrix.resize(getMaxVertexIndex());
    for (auto& v : matrix) {
        v.resize(getMaxVertexIndex(), 0);
    }
        
    for (unsigned int ei = 0; ei < edges.size(); ++ei) {
        matrix[edges[ei].first][edges[ei].second] = 1;
        matrix[edges[ei].second][edges[ei].first] = 1;
    }
    
    return matrix;
}

/*
void Dimacs::calculateGraphStats(int& maxDegree, int& minDegree, std::vector<float>& degreeHistogram) {
    auto degCopy = degrees;
    std::sort(degCopy.begin(), degCopy.end());
    size_t n = degCopy.size();
    size_t m = degreeHistogram.size();
    int cnt = 0;
    size_t di = n-1;
    for (size_t i = 0; i < m; ++i) {
        float bound = (m-i-1)*(float)n / m;
        while (degCopy[di] > bound && di > 0) {
            --di;
            ++cnt;
        }
        degreeHistogram[i] = cnt / (float)n;
        cnt = 0;
    }
    maxDegree = degCopy.front();
    minDegree = degCopy.back();
}*/

bool Dimacs::parseProblemLine(std::istringstream& ss, std::string& line) {
	if (line[0] == 'p') { // problem line
		if (!edgeNotSpecified) {
			std::string word;
			ss >> word;
			if (word != "edge") {
				error = "invalid format - edge specification missing";
				errorFlag = true;
				return false;
			}
		}
		size_t nEdges = 0;
		ss >> numVertices >> nEdges;
        // since DIMACS counts from 1 on, the storage should be larger by 1, since 0 is included implicitly
        maxVertexIndex = numVertices+1;
		edges.reserve(nEdges);
		degrees.resize(maxVertexIndex, 0);
        if (ENABLE_DEBUG_OUT_READING)
            std::cout << "problem specs: num edges = " << nEdges << ", max vertex index = " << maxVertexIndex << "\n";
	}
	return true;
}

bool Dimacs::parseSpecsLine(std::istringstream& ss, std::string& line) {
	if (line[0] == 'e') { // edge line
		unsigned int v1, v2;
		ss >> v1 >> v2;
        if (!ss) {
            // what was read in the line above does not make sense (it is very likley this is not dimacs file)
            error = "line \""+line+"\"" + "is misformed";
            errorFlag = false;
			return false;
        }
		if ((v1 > numVertices) || (v2 > numVertices) || (v1 == 0)  || (v2 == 0)) {
			error = "error detected while parsing vertices from edge spec: specs line was misformed: \""+line+"\"";
			errorFlag = true;
			return false;
		}
		//edges.push_back(std::make_pair(v1, v2));

        //adjacency list implementation
        edges[v1].push_back(v1);
        edges[v2].push_back(v2);
		degrees[v1]++;
		degrees[v2]++;
	} 
	return true;
}

/**
 * @brief 
 * @param ss
 * @param line
 * @return false on error, true otherwise 
 */
bool Dimacs::parseLine(std::istringstream& ss, std::string& line) {
    return parseProblemLine(ss, line) && parseSpecsLine(ss, line);
}

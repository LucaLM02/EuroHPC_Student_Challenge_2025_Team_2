#ifndef TEST_COMMON_HPP
#define TEST_COMMON_HPP

#include <string>
#include <sstream>
#include <iostream>

#include "common.hpp"
#include "graph.hpp"

namespace TestFunctions {
    bool IsSymmetric(Edges& edges, std::string& error_message);

    template <class Value>
    void PrintVector(std::vector<Value> to_print) {
        for ( const Value& element : to_print ) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

    template <class Value>
    std::string VecToString(std::vector<Value> vector) {
        std::stringstream ss;

        for ( const Value& element : vector ) {
            ss << element << " ";
        }

        return ss.str();
    }

    bool CheckColoring(const Graph& graph) {
        unsigned short current_color;
        std::vector<int> neighbours;
        for ( int vertex : graph.GetVertices() ) {
            current_color = graph.GetColor(vertex);
            graph.GetNeighbours(vertex, neighbours);

            if ( current_color == 0 ) {
                std::cout << "Vertex: " << vertex << " has no color" << std::endl;
                return false;
            }

            for ( int neighbour : neighbours ) {
                if ( graph.GetColor(neighbour) == current_color ) {
                    std::cout << "Vertex: " << vertex << " Color: " << current_color 
                              << " Neighbour: " << neighbour << " Color: " << graph.GetColor(neighbour) << std::endl;
                    return false;
                } 
            }
            
        }
        return true;
    }

}

#endif // TEST_COMMON_HPP

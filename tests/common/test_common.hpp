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

}

#endif // TEST_COMMON_HPP

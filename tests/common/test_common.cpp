#include "test_common.hpp"
#include <sstream>

namespace TestFunctions {
    bool IsSymmetric(Edges& edges, std::string& error_message) {
        unsigned int size=edges.size();
        std::stringstream ss;

        for (int i = 0; i < size; i++) {

            // must be 1) a matrix 2) square
            if ( edges[i].size() != size ) {
                ss << "Invalid size at row " << i << "; expected: " << size << " actual: " << edges[i].size();
                error_message = ss.str();
                return false;
            }

            // comparing elements before the diagonal element edges[i][i] 
            // with the corresponding simmetric elements
            for (int j = 0; j < i; j++) {
                if ( edges[i][j] != edges[j][i] ) {
                    ss << "Non symmetric values edges["<<i<<"]["<<j<<"]="<<edges[i][j]<<" vs edges["<<j<<"]["<<i<<"]="<<edges[j][i];
                    error_message = ss.str();
                    return false;
                }
            };
        }

        return true;
    }


}
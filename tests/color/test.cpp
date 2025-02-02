#include <iostream>

#include "color.hpp"
#include "test_common.hpp"
#include "constraints.hpp"

void test_no_constraints() {
    VertexSet vertices = { 0, 1, 3, 4, 2, 5, 7, 6}; // {1, 2, 4, 5, 3, 0, 6, 7} OK  {0, 1, 2, 3, 4, 5, 6, 7} OK
    Edges     edges(vertices.size());
    for (int i = 0; i < vertices.size(); i++) {
        edges[i].resize(vertices.size(), false);
    }
    
    edges[0][1] = true;                                                                                                     edges[0][6] = true;
    edges[1][0] = true;                     edges[1][2] = true;                     edges[1][4] = true;
                        edges[2][1] = true;                     edges[2][3] = true; edges[2][4] = true;
                                            edges[3][2] = true;                     edges[3][4] = true; edges[3][5] = true;
                        edges[4][1] = true; edges[4][2] = true; edges[4][3] = true; 
                                                                edges[5][3] = true;                                                             edges[5][7] = true;
    edges[6][0] = true;
                                                                                                        edges[7][5] = true;

    Constraints configuration(vertices.size());

    GreedyColorStrategy color_strategy;

    std::vector<unsigned short> coloring;
    unsigned short max_k;

    std::string error_message;
    if ( TestFunctions::IsSymmetric(edges, error_message) )
        std::cout << "The adjacency matrix is symmetric" << std::endl;
    else 
        std::cout << "The adjacency matrix is not symmetric" << std::endl << "error message: " << error_message << std::endl;

    color_strategy.Color(vertices, edges, configuration, 
                         coloring, max_k, vertices.size());
    
    std::cout << "Vertices: " << TestFunctions::VecToString(vertices) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(coloring) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
}

void test_with_constraints() {
    VertexSet vertices = {0, 1, 2, 3, 4, 5, 6, 7}; // { 0, 1, 3, 4, 2, 5, 7, 6} OK  {1, 2, 4, 5, 3, 0, 6, 7} OK  {0, 1, 2, 3, 4, 5, 6, 7} OK
    Edges     edges(vertices.size());
    for (int i = 0; i < vertices.size(); i++) {
        edges[i].resize(vertices.size(), false);
    }
    
    // why this shape? Well, it resembles the adjacency matrix
    edges[0][1] = true;                                                                                                     edges[0][6] = true;
    edges[1][0] = true;                     edges[1][2] = true;                     edges[1][4] = true;
                        edges[2][1] = true;                     edges[2][3] = true; edges[2][4] = true;
                                            edges[3][2] = true;                     edges[3][4] = true; edges[3][5] = true;
                        edges[4][1] = true; edges[4][2] = true; edges[4][3] = true; 
                                                                edges[5][3] = true;                                                             edges[5][7] = true;
    edges[6][0] = true;
                                                                                                        edges[7][5] = true;

    Constraints constraints(vertices.size());

    GreedyColorStrategy color_strategy;

    std::vector<unsigned short> coloring;
    unsigned short max_k;

    std::string error_message;
    if ( TestFunctions::IsSymmetric(edges, error_message) )
        std::cout << "The adjacency matrix is symmetric" << std::endl;
    else 
        std::cout << "The adjacency matrix is not symmetric" << std::endl << "error message: " << error_message << std::endl;

    color_strategy.Color(vertices, edges, constraints, 
                         coloring, max_k, vertices.size());
    
    std::cout << "Vertices: " << TestFunctions::VecToString(vertices) << std::endl;
    std::cout << "Coloring: " << TestFunctions::VecToString(coloring) << std::endl;
    std::cout << "Max k:    " << max_k << std::endl;
}

int main() {
    test_with_constraints();
    return 0;
}
#include "common.hpp"
#include "constraints.hpp"

#include <iostream>

void PrintEqualConstraints(const Constraints& c, unsigned int vertex_a) {
    std::cout << "Equal constraints of " << vertex_a << ": ";
    for (unsigned int equal_vertex_a : c.GetEqualConstraints(vertex_a)) {
        std::cout << equal_vertex_a << " ";
    }
    std::cout << std::endl;
}

void PrintDiffConstraints(const Constraints& c, unsigned int vertex_a) {
    std::cout << "Diff constraints of " << vertex_a << ": ";
    for (unsigned int equal_vertex_a : c.GetDiffConstraints(vertex_a)) {
        std::cout << equal_vertex_a << " ";
    }
    std::cout << std::endl;
}

int main() {
    VertexSet vertices = {0, 1, 2, 3, 4, 5, 6, 7};
    Constraints constraints(vertices.size());

    constraints.AddEqualConstraint(0, 1);
    std::cout << "After equal 0-1" << std::endl;
    PrintEqualConstraints(constraints, 0);
    PrintEqualConstraints(constraints, 1);

    constraints.AddEqualConstraint(0, 2);
    std::cout << "After equal 0-2" << std::endl;
    PrintEqualConstraints(constraints, 0);
    PrintEqualConstraints(constraints, 1);
    PrintEqualConstraints(constraints, 2);

    constraints.AddDiffConstraint(0, 7);
    std::cout << "After diff 0-7" << std::endl;
    PrintEqualConstraints(constraints, 0);
    PrintEqualConstraints(constraints, 1);
    PrintEqualConstraints(constraints, 2);
    PrintEqualConstraints(constraints, 7);
    // --
    PrintDiffConstraints(constraints, 0);
    PrintDiffConstraints(constraints, 1);
    PrintDiffConstraints(constraints, 2);
    PrintDiffConstraints(constraints, 7);

    constraints.AddEqualConstraint(3, 4);
    std::cout << "After equal 3-4" << std::endl;
    PrintEqualConstraints(constraints, 0);
    PrintEqualConstraints(constraints, 1);
    PrintEqualConstraints(constraints, 2);
    PrintEqualConstraints(constraints, 3);
    PrintEqualConstraints(constraints, 4);
    PrintEqualConstraints(constraints, 7);

    PrintDiffConstraints(constraints, 0);
    PrintDiffConstraints(constraints, 1);
    PrintDiffConstraints(constraints, 2);
    PrintDiffConstraints(constraints, 3);
    PrintDiffConstraints(constraints, 4);
    PrintDiffConstraints(constraints, 7);

    constraints.AddEqualConstraint(1, 3);
    std::cout << "After equal 1-3" << std::endl;
    PrintEqualConstraints(constraints, 0);
    PrintEqualConstraints(constraints, 1);
    PrintEqualConstraints(constraints, 2);
    PrintEqualConstraints(constraints, 3);
    PrintEqualConstraints(constraints, 4);
    PrintEqualConstraints(constraints, 7);
    //--
    PrintDiffConstraints(constraints, 0);
    PrintDiffConstraints(constraints, 1);
    PrintDiffConstraints(constraints, 2);
    PrintDiffConstraints(constraints, 3);
    PrintDiffConstraints(constraints, 4);
    PrintDiffConstraints(constraints, 7);
}
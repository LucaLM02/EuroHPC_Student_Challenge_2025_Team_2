#include "csr_graph.hpp"
#include "dimacs.hpp"

#include "test_common.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <chrono>


int main() {
    std::string file_name = "10_vertices_graph.clq";
    CSRGraph* graph = CSRGraph::LoadFromDimacs(file_name);
    CSRGraph* correct_graph = CSRGraph::LoadFromDimacs(file_name);

    Graph::GraphHistory history;
    history.AddAction(2, 6, Graph::GraphHistory::MERGE);
    correct_graph->MergeVertices(2,6);
    history.AddAction(2, 7, Graph::GraphHistory::ADD_EDGE);
    correct_graph->AddEdge(2,7);

    graph->AddHistory(history);

    if ( !graph->isEqual(*correct_graph) ) {
        std::cout << "Error, graphs are not equal" << std::endl;
    } else {
        std::cout << "Good job" << std::endl;
    }

    Graph::GraphHistory another_history;
    another_history.Deserialize(history.Serialize());

    CSRGraph* third_graph = CSRGraph::LoadFromDimacs(file_name);
    third_graph->AddHistory(another_history);

    if ( !third_graph->isEqual(*correct_graph) ) {
        std::cout << "Error, graphs are not equal" << std::endl;
    } else {
        std::cout << "Good job" << std::endl;
    }

}
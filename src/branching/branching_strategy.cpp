#include "branching_strategy.hpp"

RandomBranchingStrategy::RandomBranchingStrategy(int num_vertices) 
: BranchingStrategy()
{
    std::random_device dev;
    _random_generator = std::make_unique<std::mt19937>(dev());
}

std::pair<unsigned int, unsigned int> RandomBranchingStrategy::ChooseVertices(const Graph& graph) {
    do 
    {
        std::uniform_int_distribution<int> u(1, graph.GetNumVertices());
        _vertex_pair.first   = u(*_random_generator);
        _vertex_pair.second  = u(*_random_generator);

    } while ( 
        (
            graph.GetDeletedVertices().contains(_vertex_pair.first) ||
            graph.GetDeletedVertices().contains(_vertex_pair.second)
        ) &&
        graph.HasEdge(_vertex_pair.first, _vertex_pair.second) );

    return _vertex_pair;
}

DegreeBranchingStrategy::DegreeBranchingStrategy() {}

std::pair<unsigned int, unsigned int> DegreeBranchingStrategy::ChooseVertices(const Graph& graph) {
    std::pair<int, int> vertex_pair;
    const std::set<int>& deleted_vertices = graph.GetDeletedVertices();

    std::vector<int> degrees; 
    graph.GetDegrees(degrees);

    std::sort(degrees.rbegin(), degrees.rbegin());

    do 
    {
        
        

    } while (graph.HasEdge(vertex_pair.first, vertex_pair.second));
}
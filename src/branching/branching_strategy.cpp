#include "branching_strategy.hpp"

// ----------------------------- RANDOM BRANCHING STRATEGY ------------------------------

RandomBranchingStrategy::RandomBranchingStrategy(int num_vertices) 
: BranchingStrategy()
{
    std::random_device dev;
    _random_generator = std::make_unique<std::mt19937>(dev());
}


std::pair<unsigned int, unsigned int> 
RandomBranchingStrategy::ChooseVertices(const Graph &graph, PairType& type) {
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

// ----------------------------- DEGREE BRANCHING STRATEGY ------------------------------

DegreeBranchingStrategy::DegreeBranchingStrategy() {}

std::pair<unsigned int, unsigned int> 
DegreeBranchingStrategy::ChooseVertices(const Graph &graph, PairType& type) {
    std::pair<int, int> vertex_pair;
    const std::set<int>& deleted_vertices = graph.GetDeletedVertices();

    std::vector<int> degrees; 
    graph.GetDegrees(degrees);

    std::sort(degrees.rbegin(), degrees.rend());

    do 
    {
        
    } while (graph.HasEdge(vertex_pair.first, vertex_pair.second));
    return vertex_pair;
}

// ------------------------ INDEPENDENT SET BRANCHING STRATEGY --------------------------

IndependentSetBranchingStrategy::IndependentSetBranchingStrategy()
: _length_diff{1u}, _length_equal{1u}, _current_length{0u}, _is_equal_sequence{true},
  _iteration_counter{0}, _chosen_indep_set{0}
{
}

void IndependentSetBranchingStrategy::SetUpdateFrequency(unsigned int frequency)
{
    _update_frequency = frequency;
}

void IndependentSetBranchingStrategy::SetLengthDiffSequence(unsigned int length)
{
    _length_diff = length;
}

void IndependentSetBranchingStrategy::SetLengthEqualSequence(unsigned int length)
{
    _length_equal = length;
}

void IndependentSetBranchingStrategy::
    FindIndependentSets(const Graph &graph, unsigned int number)
{
    Graph& copy = *graph.Clone();
    int selected_v;

    // TODO: to implement this, maybe a Partition method on graph might be useful
    for (int counter = 0; counter < number; counter++) {
        selected_v = copy.GetVertexWithMaxDegree();
        
    }
}

std::pair<unsigned int, unsigned int> 
IndependentSetBranchingStrategy::ChooseVertices(const Graph &graph, PairType& type)
{
    // TODO: this is not enough, not all combination are explored in this way

    // 1. updates the independent sets, since some vertices have been deleted
    if ( _iteration_counter == _update_frequency ) {
        this->FindIndependentSets(graph, _independent_sets.size());
        _iteration_counter = 0;
        _chosen_indep_set = 0;
    }

    // starting again from the first set
    if ( _chosen_indep_set == _independent_sets.size() ) {
        _chosen_indep_set = 0;
    }


    std::pair<int, int> vertex_pair;
    _current_length++;
    if ( _is_equal_sequence ) {

        // TODO: così facendo però alla prima volta che itero sul set, prendo gli ultimi due e ne elimino uno.
        //       La volta successiva faccio lo stesso. Però uno degli ultimi 2 è quello rimasto dalla volta precedente
        //       Andrebbe modificato
        //       Ancora peggio nel caso in cui faccio `diff`

        // popping another element from the same set (popping since it will be merged with the first one, disappearing)
        vertex_pair.second = _independent_sets[_chosen_indep_set][_independent_sets[_chosen_indep_set].size()-1];
        _independent_sets[_chosen_indep_set].pop_back();
        
        vertex_pair.first = _independent_sets[_chosen_indep_set][_independent_sets[_chosen_indep_set].size()-1];

        _chosen_indep_set++;

        if ( _current_length == _length_equal ) {
            _current_length = 0;
            _is_equal_sequence = false;
        }
    } else {

        bool edge_found = false;

        int i = _independent_sets[_chosen_indep_set].size() - 1;
        int second_set = _chosen_indep_set + 1;

        while ( !edge_found ) {
            // if no edge can be created from chosen set to the other, changing second set
            if ( i == 0 ) {
                i = _independent_sets[_chosen_indep_set].size() - 1;
                second_set++;
            }

            // choosing second set
            if ( second_set == _independent_sets.size() ) {
                second_set = 0;
            }

            for (int j = 0; j < _independent_sets[second_set].size() - 1; j++) {
                // getting element from first chosen set
                vertex_pair.first = _independent_sets[_chosen_indep_set][i];

                // getting element from second chosen set
                vertex_pair.second = _independent_sets[second_set][j];
                
                // quitting if a valid edge is found
                if ( !graph.HasEdge(vertex_pair.first, vertex_pair.second) ) {
                    edge_found = true;
                    break;
                }
            }

            i--;

        }

        _chosen_indep_set++;

        if ( _current_length == _length_diff ) {
            _current_length = 0;
            _is_equal_sequence = true;
        }
    }

    _iteration_counter++;
    return vertex_pair;
}

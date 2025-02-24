#include "branching_strategy.hpp"
#include <algorithm>

// ----------------------------- RANDOM BRANCHING STRATEGY ------------------------------

RandomBranchingStrategy::RandomBranchingStrategy() 
: BranchingStrategy()
{
    std::random_device dev;
    _random_generator = std::make_unique<std::mt19937>(dev());
}


std::pair<int, int> 
RandomBranchingStrategy::ChooseVertices(Graph &graph) {
    //check if the graph is complete
    int n = graph.GetNumVertices();
    if(graph.GetNumEdges() == n*(n-1)/2) {
        return std::make_pair(-1, -1);
    }
    do 
    {
        std::uniform_int_distribution<int> u(0, n-1);
        _vertex_pair.first   = graph.GetVertices()[u(*_random_generator)];
        _vertex_pair.second  = graph.GetVertices()[u(*_random_generator)];

    } while ( 
        /*
        (
            graph.GetDeletedVertices().contains(_vertex_pair.first) ||
            graph.GetDeletedVertices().contains(_vertex_pair.second)
        ) &&
         */
        graph.HasEdge(_vertex_pair.first, _vertex_pair.second) || _vertex_pair.first == _vertex_pair.second );

    return _vertex_pair;
}

// ----------------------------- DEGREE BRANCHING STRATEGY ------------------------------

DegreeBranchingStrategy::DegreeBranchingStrategy() {}

std::pair<int, int> 
DegreeBranchingStrategy::ChooseVertices(Graph &graph) {
    std::pair<int, int> vertex_pair;

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
/* void IndependentSetBranchingStrategy::FindIndependentSets(Graph &graph, unsigned int number)
{
    _independent_sets.clear();
    std::unordered_set<int> remaining_vertices(graph.GetVertices().begin(), graph.GetVertices().end());

    while (!remaining_vertices.empty()) {
        std::vector<int> independent_set;
        for (auto it = remaining_vertices.begin(); it != remaining_vertices.end(); ) {
            bool is_independent = true;
            for (int v : independent_set) {
                if (graph.HasEdge(v, *it)) {
                    is_independent = false;
                    break;
                }
            }
            if (is_independent) {
                independent_set.push_back(*it);
                it = remaining_vertices.erase(it);
            } else {
                ++it;
            }
        }
        _independent_sets.push_back(independent_set);
    }
}
*/

std::pair<int, int> 
IndependentSetBranchingStrategy::ChooseVertices(Graph &graph)
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

/* std::pair<int, int> IndependentSetBranchingStrategy::ChooseVertices(Graph &graph)
{
    if (_iteration_counter == _update_frequency || _independent_sets.empty()) {
        this->FindIndependentSets(graph, _independent_sets.size());
        _iteration_counter = 0;
        _chosen_indep_set = 0;
    }

    if (_chosen_indep_set >= _independent_sets.size()) {
        _chosen_indep_set = 0;
    }

    std::pair<int, int> vertex_pair;
    _current_length++;
    if (_is_equal_sequence) {
        if (_independent_sets[_chosen_indep_set].size() < 2) {
            _chosen_indep_set++;
            return ChooseVertices(graph);
        }
        vertex_pair.second = _independent_sets[_chosen_indep_set].back();
        _independent_sets[_chosen_indep_set].pop_back();
        vertex_pair.first = _independent_sets[_chosen_indep_set].back();
        _chosen_indep_set++;

        if (_current_length == _length_equal) {
            _current_length = 0;
            _is_equal_sequence = false;
        }
    } else {
        for (size_t i = 0; i < _independent_sets.size(); i++) {
            for (size_t j = i + 1; j < _independent_sets.size(); j++) {
                for (int v1 : _independent_sets[i]) {
                    for (int v2 : _independent_sets[j]) {
                        if (!graph.HasEdge(v1, v2)) {
                            vertex_pair = {v1, v2};
                            _chosen_indep_set++;
                            _iteration_counter++;
                            return vertex_pair;
                        }
                    }
                }
            }
        }
    }

    _iteration_counter++;
    return vertex_pair;
}
*/

std::pair<int, int> NeighboursBranchingStrategy::ChooseVertices(Graph &graph)
{
    std::vector<int> vertices = graph.GetVertices();
    int vertex_x = -1, vertex_y = -1, vertex_w, vertex_z;
    std::set<int> w_neighbours;
    std::set<int> z_neighbours;

    int max_common_neighbours = -1;
    int curr_common_neighbours;

    for ( int i = 0; i < vertices.size(); i++ ) {
        vertex_w = vertices[i];
        graph.GetNeighbours(vertex_w, w_neighbours);

        for ( int j = i + 1; j < vertices.size(); j++ ) {
            vertex_z = vertices[j];
            
            // skipping adjacent vertices
            if ( graph.HasEdge(vertex_w, vertex_z) ) {
                continue;
            }

            graph.GetNeighbours(vertex_z, z_neighbours);

            curr_common_neighbours = 0;
            for ( int neighbour_z : z_neighbours ) {
                if ( w_neighbours.contains(neighbour_z) ) {
                    curr_common_neighbours++;
                }
            }

            if ( max_common_neighbours < curr_common_neighbours ) {
                max_common_neighbours = curr_common_neighbours;
                vertex_x = vertex_w;
                vertex_y = vertex_z;
            }
        }
    }
    

    return {vertex_x, vertex_y};
}

/* std::pair<int, int> NeighboursBranchingStrategy::ChooseVertices(Graph &graph)
{
    std::vector<int> vertices = graph.GetVertices();
    int vertex_x = -1, vertex_y = -1;
    int max_common_neighbours = -1;

    for (size_t i = 0; i < vertices.size(); i++) {
        for (size_t j = i + 1; j < vertices.size(); j++) {
            if (graph.HasEdge(vertices[i], vertices[j])) continue;

            int common_neighbours = graph.CountCommonNeighbours(vertices[i], vertices[j]); // consider implementing CountCommonNeighbours
            if (common_neighbours > max_common_neighbours) {
                max_common_neighbours = common_neighbours;
                vertex_x = vertices[i];
                vertex_y = vertices[j];
            }
        }
    }

    return {vertex_x, vertex_y};
}
*/


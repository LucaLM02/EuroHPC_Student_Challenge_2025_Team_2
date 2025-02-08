#include "dsatur_color.hpp"

void DSaturColorStrategy::Color(Graph &graph, unsigned short &max_k) const
{
    std::vector<unsigned short> coloring(graph.GetHighestVertex() + 1);

    DSaturList list(graph);
    int selected_vertex;
    unsigned short selected_color;

    std::vector<int> neighbours;

    max_k = 0;
    while ( !list.IsEmpty() ) {
        // retrieves the element with highest degree among the ones with highest 
        // saturation degree
        selected_vertex = list.PopHighestVertex();

        // colors the vertex in a greedy fashion
        selected_color  = GreedyFindColor(graph, selected_vertex, coloring, max_k + 1);
        coloring[selected_vertex] = selected_color;

        // updates the maximum color
        if ( selected_color > max_k ) {
            max_k = selected_color;
        }

        // updating the saturation degree list
        graph.GetNeighbours(selected_vertex, neighbours);
        for ( int neighbour : neighbours ) {
            if ( coloring[neighbour] > 0 ) {
                continue;
            }
            list.AddNeighbourColor(neighbour, selected_color);
        }
    }

    graph.SetFullColoring(coloring);
    
}


// ================================= DSATURLIST PUBLIC ===================================


DSaturList::DSaturList(const Graph& graph)
: _vertex_to_item(graph.GetHighestVertex()+1),   // + 1 because highest vertex must be included
  _vertex_to_neighbour_colors(graph.GetHighestVertex()+1),
  _sat_degree_to_list(1, nullptr)
{
    std::vector<int> vertices = graph.GetVertices();
    std::vector<int> degrees = graph.GetFullDegrees();

    static auto ascendingCompare = 
    [&](int v, int w) -> bool {
        return degrees[v] < degrees[w];
    };

    std::sort(vertices.begin(), vertices.end(), ascendingCompare);

    DSaturItem* prev = nullptr;
    DSaturItem* item = new DSaturItem(0, degrees[vertices[0]], vertices[0], nullptr);
    _sat_degree_to_list[0] = item;          // all saturation degrees set to 0
    _vertex_to_item[vertices[0]] = item;    // saving first item
    prev = item;

    // removing first vertex
    vertices.erase(vertices.begin());

    // looping over the remaining vertices
    for (int vertex : vertices) {
        item = new DSaturItem(0, degrees[vertex], vertex, prev);
        prev->next = item;
        prev = item;
        _vertex_to_item[vertex] = item;
    }

    item->next = nullptr;
    _last_item = item;
    _last_degree = item->sat_degree;
    
}

void DSaturList::AddNeighbourColor(int vertex, unsigned short color)
{
    if ( _vertex_to_neighbour_colors[vertex].insert(color).second ) {
        IncreaseSatDegree(vertex);
    }
}

int DSaturList::GetLowestSatDegree() const
{
    if ( this->IsEmpty() ) {
        return -1;
    }
    int vertex = GetLowestVertex();
    return _vertex_to_item[vertex]->sat_degree;
}

int DSaturList::GetLowestVertex() const
{
    if ( this->IsEmpty() ) {
        return -1;
    }

    int first_degree_counter;
    
    // searching the first item
    for (first_degree_counter = 0; 
         first_degree_counter <= _last_degree &&
         _sat_degree_to_list[first_degree_counter] == nullptr;
         first_degree_counter++) {}

    if ( _sat_degree_to_list[first_degree_counter] == nullptr ) {
        return -1;
    }
    return _sat_degree_to_list[first_degree_counter]->vertex;
}

int DSaturList::PopLowestVertex()
{
    int first_degree_counter;

    // searching the first item
    for (first_degree_counter = 0; 
         _sat_degree_to_list[first_degree_counter] == nullptr;
         first_degree_counter++) {}

    // removing the first item
    DSaturItem* first_item = _sat_degree_to_list[first_degree_counter];
    _sat_degree_to_list[first_degree_counter] = first_item->next;
    if ( first_item->next != nullptr ) {
        first_item->next->prev = nullptr;
    }

    int removed_vertex = first_item->vertex;
    free(first_item);

    return removed_vertex;
}

int DSaturList::GetHighestSatDegree() const
{
    return _last_degree;
}

int DSaturList::GetHighestVertex() const
{
    if ( _last_item == nullptr ) {
        return -1;
    }
    return _last_item->vertex;
}

int DSaturList::PopHighestVertex()
{
    int removed_vertex = _last_item->vertex;
    DSaturItem* removed_item = _last_item;

    if ( _last_item->prev != nullptr ) {
        _last_item = _last_item->prev;
        _last_item->next = nullptr;      // it cannot have another element after
    } else {
        // in this case it is the first of a list, searching the new _last_item in the
        // previous lists

        _sat_degree_to_list[_last_degree] = nullptr;

        // scrolling all the lists
        int last_degree_counter = _last_degree;
        for (last_degree_counter = _last_degree; 
            last_degree_counter >= 0 &&
            _sat_degree_to_list[last_degree_counter] == nullptr;
            last_degree_counter--) {}
        
        if ( last_degree_counter < 0 ) {
            // list has been cleared
            _last_item = nullptr;
            _last_degree = -1;
            free(removed_item);
            return removed_vertex;
        }

        _last_degree = last_degree_counter;

        // finding the last element of the list
        _last_item = _sat_degree_to_list[last_degree_counter];

        if ( _last_item == nullptr ) {
            return true;
        }

        while ( _last_item->next != nullptr ) {
            _last_item = _last_item->next;
        }
    }

    free(removed_item);
    return removed_vertex;
}

bool DSaturList::IsEmpty() const
{
    return _last_item == nullptr;
}

const DSaturItem *DSaturList::operator[](int degree) const
{
    return _sat_degree_to_list[degree];
}

// ================================ DSATURLIST PROTECTED =================================

bool DSaturList::IncreaseSatDegree(int vertex, int increment)
{
    DSaturItem* item = _vertex_to_item[vertex];

    int prev_sat_degree = item->sat_degree;
    int new_sat_degree = item->sat_degree + increment;

    if ( new_sat_degree < 0 ) {
        return false;
    }

    item->sat_degree = new_sat_degree;

    // removing from current list
    if ( item->next != nullptr ) {
        item->next->prev = item->prev;
    }

    if ( item->prev != nullptr ) {
        item->prev->next = item->next;
    } else {
        // moving the first element
        _sat_degree_to_list[prev_sat_degree] = item->next;
    }


    // creating new saturation degree lists, if needed
    if ( new_sat_degree >= _sat_degree_to_list.size() ) {
        _sat_degree_to_list.resize(new_sat_degree + 1, nullptr);
    }

    if ( new_sat_degree > _last_degree ) {
        // updating the _last_item in case the new highest saturation degree list has
        // been created
        _last_item = item;
        _last_degree = new_sat_degree;

        // this imply that item is easy to insert
        _sat_degree_to_list[new_sat_degree] = item;
        item->next = nullptr;
        item->prev = nullptr;
    } else if ( new_sat_degree == _last_degree && 
                item->degree > _last_item->degree ) {
        // N.B.: using > and not >= because if item X with degree A is present and a 
        //       new item Y with degree A is added, then Y will come first than X
        //       So it will never be a new _last_item
        // updating the _last_item in case `item` it is not in a branch new saturation 
        // degree list but it is appended at the end of the highest saturation degree
        // list
        _last_item->next = item;
        item->prev = _last_item;
        item->next = nullptr;
        _last_item = item;
    } else {
        // in this case the new item is not a new _last_item

        DSaturItem* new_prev = _sat_degree_to_list[new_sat_degree];
        if ( new_prev == nullptr ) {
            // if list is empty
            _sat_degree_to_list[new_sat_degree] = item;
            item->prev = nullptr;
            item->next = nullptr;
        } else if ( new_prev->degree >= item->degree ) {
            // if list is non-empty and item is placed as first
            _sat_degree_to_list[new_sat_degree] = item;
            new_prev->prev = item;
            item->next = new_prev;
            item->prev = nullptr;
        } else {
            // if list is non-empty and item is not placed as first then scrolling until
            // either the end or an item with higher degree is found
            while ( new_prev->next != nullptr && new_prev->next->degree < item->degree ) {
                new_prev = new_prev->next;
            }

            // if end is reached, handling nullptr
            if ( new_prev->next == nullptr ) {
                new_prev->next = item;
                item->next = nullptr;
                item->prev = new_prev;
            } else {
                // otherwise, modifying both new_prev and new_next items
                new_prev->next->prev = item;
                item->next = new_prev->next;
                new_prev->next = item;
                item->prev = new_prev;
            }
            
        }
    }


    return true;
}

bool DSaturList::DecreaseSatDegree(int vertex, int decrement)
{
    return this->IncreaseSatDegree(vertex, -decrement);
}

bool DSaturList::IncreaseDegree(int vertex, unsigned int increment)
{
    DSaturItem* item = _vertex_to_item[vertex];
    int new_degree = item->degree + increment;

    // checking if the order remains the same
    if ( item->next == nullptr || item->next->degree >= new_degree ) {
        return true;
    }

    // removing the item
    if ( item->prev == nullptr ) {
        _sat_degree_to_list[item->sat_degree] = item->next;
        item->next->prev = nullptr;
    } else {
        item->prev->next = item->next;
        item->next->prev = item->prev;
    }

    // scrolling the list
    DSaturItem* new_prev = item;
    while ( new_prev->next != nullptr && new_prev->next->degree < new_degree ) {
        new_prev = new_prev->next;
    }

    // inserting the item
    item->prev = new_prev;
    if ( new_prev->next == nullptr ) {
        new_prev->next = item;
        item->next = nullptr;
    } else {
        item->next = new_prev->next;
        new_prev->next = item;
        item->next->prev = item;
    }

    return true;
}

bool DSaturList::DecreaseDegree(int vertex, unsigned int decrement)
{

    DSaturItem* item = _vertex_to_item[vertex];
    int new_degree = item->degree - decrement;

    if ( new_degree < 0 ) {
        return false;
    }

    // checking if the order remains the same
    if ( item->prev == nullptr || item->prev->degree <= new_degree ) {
        return true;
    }

    // removing the item
    if ( item->next == nullptr ) {
        item->prev->next = nullptr;
    } else {
        item->prev->next = item->next;
        item->next->prev = item->prev;
    }

    // scrolling the list
    DSaturItem* new_next = item;
    while ( new_next->prev != nullptr && new_next->prev->degree > new_degree ) {
        new_next = new_next->prev;
    }

    // inserting the item
    item->next = new_next;
    if ( new_next->prev == nullptr ) {
        new_next->prev = item;
        _sat_degree_to_list[item->sat_degree] = item;
    } else {
        item->prev = new_next->prev;
        new_next->prev = item;
        item->prev->next = item;
    }

    return true;
}

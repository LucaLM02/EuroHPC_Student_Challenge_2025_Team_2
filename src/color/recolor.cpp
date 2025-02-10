#include "recolor.hpp"

unsigned int GreedySwapRecolorStrategy::Recolor(Graph& graph) const 
{
    std::vector<unsigned short> coloring = graph.GetFullColoring();
    SwapRecolorStructure data(graph, coloring, 50);
    data.FillWithData();
    int result = data.Recolor();

    graph.SetFullColoring(coloring);
    return result;
}

VertexRecolorData::VertexRecolorData()
    : _old_color{NOT_ASSIGNED},
      _vertex{NOT_ASSIGNED}
{
    std::random_device dev;
    _random_generator = std::make_unique<std::mt19937>(dev());
}

VertexRecolorData::VertexRecolorData(int vertex,
                                     std::vector<unsigned short>* coloring,
                                     unsigned short max_color)
    : _vertex{vertex}, _old_color{NOT_ASSIGNED},
      _coloring{coloring}, _max_color{max_color}
{
    std::random_device dev;
    _random_generator = std::make_unique<std::mt19937>(dev());
}

inline void VertexRecolorData::AssignColor(unsigned short color)
{
    if ( _old_color == NOT_ASSIGNED ) {
        _old_color = (*_coloring)[_vertex];
    }
    (*_coloring)[_vertex] = color;
}

SwapRecolorStructure::SwapRecolorStructure(Graph &graph, std::vector<unsigned short>& coloring, int threshold)
: _graph{graph}, _vertex_to_data(_graph.GetHighestVertex()+1), 
  _dont_color{false}, _threshold{threshold}, _coloring{coloring}
{
}

void SwapRecolorStructure::FillWithData()
{
    _graph.SortByColor();
    std::vector<int> vertices = _graph.GetVertices();
    unsigned int max_color = _graph.GetColor(vertices[0]);


    int counter = 0;
    for (int vertex : vertices ) {
        if ( _graph.GetColor(vertex) == max_color ) {
            counter++;
            if ( counter == _threshold ) {
                _dont_color = true;
                return;
            }
        } else {
            break;
        }
    }

    std::vector<int> neighbours;
    std::vector<int> other_neighbours;
    
    /**
     * @brief for each colors, true if a particular vertix can be recolored with it
     *        false otherwise
     * @note + 1 since first element has index 0 and is skipped
     */
    std::vector<bool> colors_availability(max_color + 1);   
    /**
     * @brief all colors with colors_availability[color] = true
     */
    std::vector<unsigned short> available_colors;
    available_colors.reserve(max_color);

    std::vector<unsigned short> full_coloring = _graph.GetFullColoring();

    for ( int vertex : vertices ) {
        if ( full_coloring[vertex] < max_color ) {
            break;
        }

        // initilizing in _vertex_to_data each neighbour of the current vertex
        _graph.GetNeighbours(vertex, neighbours);
        // also initilizing the current vertex data structure
        //neighbours.push_back(vertex);
        for ( int neighbour : neighbours ) {

            VertexRecolorData& data = _vertex_to_data[neighbour];
            // if is assigned then this vertex was already visited and added
            if ( data.GetVertex() == VertexRecolorData::NOT_ASSIGNED ) 
            {
                _graph.GetNeighbours(neighbour, other_neighbours);
                data.InitColoring(&_coloring, max_color);
                data.InitVertex(neighbour, other_neighbours);
            }
        }
    }
}

bool SwapRecolorStructure::Recolor()
{
    if ( _dont_color ) {
        return false;
    }

    std::vector<int> vertices = _graph.GetVertices();
    unsigned int max_color = _graph.GetColor(vertices[0]);

    // finds the first vertex with color < max_color
    int index;
    for (index = 0; 
         index < vertices.size() && _graph.GetColor(vertices[index]) == max_color; 
         index++)
    {
    }

    // erases all vertices with non maximum color
    vertices.erase(vertices.begin()+index, vertices.end());

    return this->RecolorBody(vertices);
}

bool SwapRecolorStructure::RecolorBody(std::vector<int> &vertices)
{
    // recursion basic case
    if ( vertices.size() == 0 ) {
        return true;
    }

    unsigned short max_color = _graph.GetColor(vertices[0]);
    // removes the current vertex
    int current_vertex = vertices[vertices.size()-1];
    vertices.pop_back();

    std::vector<int> neighbours;
    _graph.GetNeighbours(current_vertex, neighbours);
    std::vector<bool> has_been_recolored(neighbours.size());

    bool successfully_recolored;

    VertexRecolorData* current_neighbour_data;

    /* Handled cases:
        1) v1 <----> neighb <----> v2
        
        2)     v1              v2
               |               |
           neighb_1         neighb_2

        3)     v1              v2
               |               |
           neighb_1 <----> neighb_2            
    */

    // for each color, tries to assign it to the current vertex and change the neighbours with that color.
    // The neighbours with a different color cannot be colored to the current color until the current vertex
    // changes it's color.
    // Then, when a valid coloring is found, starts coloring the next vertex (recursive call)
    for (unsigned short color=1; color < max_color; color++ ) {

        _coloring[current_vertex] = color;
        successfully_recolored = true;

        std::fill(has_been_recolored.begin(), has_been_recolored.end(), false);

        for (int neighbour_index = 0; neighbour_index < neighbours.size(); neighbour_index++ ) {
            current_neighbour_data = &(_vertex_to_data[neighbours[neighbour_index]]);


            if ( current_neighbour_data->GetCurrentColor() == color ) {
                // if it has the same color, then change it

                has_been_recolored[neighbour_index] = true;

                // if unable to recolor, reverts the actions
                if ( !current_neighbour_data->Recolor() ) {
                    for ( int i=0; i < neighbour_index; i++ ) {
                        if ( has_been_recolored[i] ) {
                            // reverts color to what was previously
                            _vertex_to_data[neighbours[i]].AssignColor(color);
                        } 
                    }
                    // going to the next color
                    successfully_recolored = false;
                    break;
                }
            } 
        }

        // going to the next color
        if ( !successfully_recolored ) {
            _coloring[current_vertex] = max_color;
            continue;
        }

        // recolors also the next vertices
        if ( this->RecolorBody(vertices) ) {
            vertices.push_back(current_vertex);
            return true;
        } else {
            // reverting all the actions done during this for-each color iteration

            for (int neighbour_index = 0; neighbour_index < neighbours.size(); neighbour_index++ ) {
                if ( has_been_recolored[neighbour_index] ) {
                    // reverts color to what was previously
                    _vertex_to_data[neighbours[neighbour_index]].AssignColor(color);
                } 
            }
        }
    }

    vertices.push_back(current_vertex);
    return false;
}

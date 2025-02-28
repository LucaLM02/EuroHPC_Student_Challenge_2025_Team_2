#include "graph.hpp"

void Graph::GraphHistory::AddAction(int v, int u, bool action)
{
    _vertices.push_back(std::pair(v, u));
    _actions.push_back(action);
}

std::string Graph::GraphHistory::Serialize() const
{
    std::ostringstream oss;
    oss << _vertices.size() << " ";
    for ( int i=0; i < _vertices.size(); i++ ) {
        oss << _vertices[i].first << " " << _vertices[i].second << " " << _actions[i] << " ";
    }
    return oss.str();
}

void Graph::GraphHistory::Deserialize(const std::string &data)
{
    std::istringstream iss(data);
    size_t size;
    
    iss >> size;
    _vertices.resize(size);
    _actions.resize(size);
    for ( int i=0; i < size; i++ ) {
        iss >> _vertices[i].first;
        iss >> _vertices[i].second;

        int action;
        iss >> action;
        _actions[i] = (action == 1);
    }
}

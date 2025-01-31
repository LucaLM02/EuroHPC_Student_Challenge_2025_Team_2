#include "constraints.hpp"

Constraints& Constraints::operator=(const Constraints& other) {
    _equal_constraints = other._equal_constraints;
    _diff_constraints = other._diff_constraints;

    return *this;
}

/* -------------------- EQUAL CONSTRAINTS ---------------------*/
const std::set<unsigned int>& Constraints::GetEqualConstraints(unsigned int vertex_a) const{
    return _equal_constraints.at(vertex_a);
}
void Constraints::AddEqualConstraint(unsigned int vertex_a, unsigned int vertex_b) {

    // buffers needed since the sets _*_constraints[vertex_a/b] are modified during the
    // for-each on them
    std::set<unsigned int> _equal_constraints_a = _equal_constraints[vertex_a];
    std::set<unsigned int> _equal_constraints_b = _equal_constraints[vertex_b];
    std::set<unsigned int> _diff_constraints_a  =_diff_constraints[vertex_a];
    std::set<unsigned int> _diff_constraints_b  =_diff_constraints[vertex_b];

    for (unsigned int equal_vertex_a : _equal_constraints_a) {
        // [equal rule 1]
        _AtomicAddEqualConstraint(equal_vertex_a, vertex_b);

        for (unsigned int equal_vertex_b : _equal_constraints_b) {
            // [equal rule 2]
            _AtomicAddEqualConstraint(equal_vertex_a, equal_vertex_b);
        }

        for (unsigned int diff_vertex_b : _diff_constraints_b) {
            // [diff rule 2]
            _AtomicAddDiffConstraint(equal_vertex_a, diff_vertex_b);
        }
    }

    for (unsigned int equal_vertex_b : _equal_constraints_b) {
        // [equal rule 1]
        _AtomicAddEqualConstraint(vertex_a, equal_vertex_b);

        for (unsigned int diff_vertex_a : _diff_constraints_a) {
            // [diff rule 2]
            _AtomicAddDiffConstraint(diff_vertex_a, equal_vertex_b);
        }
    }

    for (unsigned int diff_vertex_a : _diff_constraints_a) {
        // [diff rule 1]
        _AtomicAddDiffConstraint(diff_vertex_a, vertex_b);
    }

    for (unsigned int diff_vertex_b : _diff_constraints_b) {
        // [diff rule 1]
        _AtomicAddDiffConstraint(vertex_a, diff_vertex_b);
    }

    _AtomicAddEqualConstraint(vertex_a, vertex_b);
}

void Constraints::AddEqualConstraint(std::pair<unsigned int, unsigned int> vertices) {
    AddEqualConstraint(vertices.first, vertices.second);
}

void Constraints::RemoveEqualConstraint(unsigned int vertex_a, unsigned int vertex_b) {
    _equal_constraints[vertex_a].erase(vertex_b);
    _equal_constraints[vertex_b].erase(vertex_a);
}

void Constraints::RemoveEqualConstraint(std::pair<unsigned int, 
                                                  unsigned int> vertices) {
    _equal_constraints[vertices.first].erase(vertices.second);
    _equal_constraints[vertices.second].erase(vertices.first);
}

void Constraints::RemoveEqualConstraints(unsigned int vertex_a) {
    _equal_constraints[vertex_a].clear();
}

bool Constraints::IsEqualConstraint(unsigned int vertex_a, 
                                    unsigned vertex_b) const {
    return _equal_constraints.at(vertex_a).contains(vertex_b);
}

bool Constraints::IsEqualConstraint(std::pair<unsigned int, 
                                              unsigned int> vertices) const {
    return _equal_constraints.at(vertices.first).contains(vertices.second);
}

/* -------------------- DIFF CONSTRAINTS ----------------------*/
const std::set<unsigned int>& Constraints::GetDiffConstraints(
                                                unsigned int vertex_a) const {
    return _diff_constraints.at(vertex_a);
}

void Constraints::AddDiffConstraint(unsigned int vertex_a, unsigned int vertex_b) {
    _AtomicAddDiffConstraint(vertex_a, vertex_b);

    // no buffers needed for _*_constraints_a/b because they are not modified
    // during the for loop

    for (unsigned int equal_vertex_a : _equal_constraints[vertex_a]) {
        // [equal rule 1]
        _AtomicAddDiffConstraint(equal_vertex_a, vertex_b);

        for (unsigned int equal_vertex_b : _equal_constraints[vertex_b]) {
            // [equal rule 2]
            _AtomicAddDiffConstraint(equal_vertex_a, equal_vertex_b);
        }
    }

    for (unsigned int equal_vertex_b : _equal_constraints[vertex_b]) {
        // [equal rule 1]
        _AtomicAddDiffConstraint(vertex_a, equal_vertex_b);

    }

}

void Constraints::AddDiffConstraint(std::pair<unsigned int, 
                                              unsigned int> vertices) {
    _diff_constraints[vertices.first].insert(vertices.second);
    _diff_constraints[vertices.second].insert(vertices.first);
}

void Constraints::RemoveDiffConstraint(unsigned int vertex_a, unsigned int vertex_b) {
    _diff_constraints[vertex_a].erase(vertex_b);
    _diff_constraints[vertex_b].erase(vertex_a);
}

void Constraints::RemoveDiffConstraint(std::pair<unsigned int, 
                                                 unsigned int> vertices) {
    _diff_constraints[vertices.first].erase(vertices.second);
    _diff_constraints[vertices.second].erase(vertices.first);
}

void Constraints::RemoveDiffConstraints(unsigned int vertex_a) {
    _diff_constraints[vertex_a].clear();
}

bool Constraints::IsDiffConstraint(unsigned int vertex_a, 
                                   unsigned vertex_b) const {
    return _diff_constraints.at(vertex_a).contains(vertex_b);
}

bool Constraints::IsDiffConstraint(std::pair<unsigned int, 
                                            unsigned int> vertices) const {
    return _diff_constraints.at(vertices.first).contains(vertices.second);
}

void Constraints::_AtomicAddEqualConstraint(unsigned int vertex_a, 
                                            unsigned int vertex_b) {
    _equal_constraints[vertex_a].insert(vertex_b);
    _equal_constraints[vertex_b].insert(vertex_a);
}
void Constraints::_AtomicAddDiffConstraint(unsigned int vertex_a, 
                                           unsigned int vertex_b) {
    _diff_constraints[vertex_a].insert(vertex_b);
    _diff_constraints[vertex_b].insert(vertex_a);
}

/* -------------------- OTHER UTILITIES -----------------------*/
void Constraints::Swap() {
    VertexMap<std::set<unsigned int>> tmp(_equal_constraints);
    _equal_constraints = _diff_constraints;
    _diff_constraints = tmp;
}
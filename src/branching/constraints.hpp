#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include "common.hpp"
#include <map>
#include <set>

/*
    IMPORTANT NOTES:
        - isn't it better to include neighbous in the list of diff-constraints? In this way
          I have to cycle only through one list and when I add an equal-constraint, the other
          vertex will inherit my neighbours as diff-constraints.
          More memory consumption but faster because that vertex has all the constraints to
          respect "at a glance"
        - the code is not correct yet. Need to implement the merge of the constraitns when a new
          constraint is added
        - 

*/

/*
    @brief container of the constraints which are applied to vertices colors

    @details
    Given 2 vertices vertex_a and vertex_b, constraints can be of 2 types:
        - vertex_a and vertex_b must have the same color (*EqualConstraint)
        - vertex_a and vertex_b must have different colors (*DiffConstraint)

    Note that adding a new (equal- or diff-)constraint will trigger an update of constraints of
    the 2 vertices involved and all their (equal- or diff-)constrained vertices.
    For more details look at the documentations of AddEqualConstraint and AddDiffConstraint

    Simply, constraints can be added (Add*), removed (Remove*) or looked up (Is*Constraint). 
    It is possible to access all the constraints of a certain type that a vertex has (Get*Constraints).
*/
class Constraints {
    public:
        /*
            @brief empty constructor, used usually if the variable of this type is then filled with the assignment-operator
        */
        Constraints() 
        {

        }

        Constraints(unsigned int num_vertices) 
        {
            for (int i = 0; i < num_vertices; i++) {
                _equal_constraints[i] = std::set<unsigned int>();
                _diff_constraints[i]  = std::set<unsigned int>();
            }
        }

        Constraints(Constraints &other) 
        : _equal_constraints(other._equal_constraints), 
          _diff_constraints(other._diff_constraints)
        {

        }

        Constraints& operator=(const Constraints& other);

        /* -------------------- EQUAL CONSTRAINTS ---------------------*/
        const std::set<unsigned int>& GetEqualConstraints(unsigned int vertex_a) const;

        /*
            @brief add an equal constraint between vertex_a and vertex_b. Updates accordingly the other constraints
            @details 
            First of all, it adds an equal-constraint between vertex_a and vertex_b.
            Then, it adds other constraints according to the following rules:
                [diff rule 1 ] color_a != color_c => color_b != color_c     
                [diff rule 2 ] color_a != color_c and color_b = color_d => color_d != color_c 
                [equal rule 1] color_a  = color_c => color_b = color_c
                [equal rule 2] color_a  = color_c and color_b = color_d => color_c = color_d
            They are also valid if a and b are swapped.

            Note that AddEqualConstraint(a, b); RemoveEqualConstraint(a, b) does not (necessarily) 
            revert the action, since the constrains generated following the above rules cannot be 
            removed since they are note tracked.
            A child class of this class could implement this.

            @param vertex_a first vertex constrained
            @param vertex_b second vertex constrained
        */
        void AddEqualConstraint(unsigned int vertex_a, unsigned int vertex_b);

        void AddEqualConstraint(std::pair<unsigned int, unsigned int> vertices);

        void RemoveEqualConstraint(unsigned int vertex_a, unsigned int vertex_b);

        void RemoveEqualConstraint(std::pair<unsigned int, unsigned int> vertices);

        void RemoveEqualConstraints(unsigned int vertex_a);

        bool IsEqualConstraint(unsigned int vertex_a, unsigned vertex_b) const;

        bool IsEqualConstraint(std::pair<unsigned int, unsigned int> vertices) const;

        /* -------------------- DIFF CONSTRAINTS ----------------------*/
        const std::set<unsigned int>& GetDiffConstraints(unsigned int vertex_a) const;

        /*
            @brief add a diff constraint between vertex_a and vertex_b. Updates accordingly the other constraints
            @details 
            First of all, it adds a diff-constraint between vertex_a and vertex_b.
            Then, it adds other constraints according to the following rules:
                [equal rule 1] color_a  = color_c => color_b != color_c
                [equal rule 2] color_a  = color_c and color_b = color_d => color_c != color_d
            They are also valid if a and b are swapped.

            Note that AddDiffConstraint(a, b); RemoveDiffConstraint(a, b) does not (necessarily) 
            revert the action, since the constrains generated following the above rules cannot be 
            removed since they are note tracked.
            A child class of this class could implement this.

            @param vertex_a first vertex constrained
            @param vertex_b second vertex constrained
        */
        void AddDiffConstraint(unsigned int vertex_a, unsigned int vertex_b);

        void AddDiffConstraint(std::pair<unsigned int, unsigned int> vertices);

        void RemoveDiffConstraint(unsigned int vertex_a, unsigned int vertex_b);

        void RemoveDiffConstraint(std::pair<unsigned int, unsigned int> vertices);

        void RemoveDiffConstraints(unsigned int vertex_a);

        bool IsDiffConstraint(unsigned int vertex_a, unsigned vertex_b) const;

        bool IsDiffConstraint(std::pair<unsigned int, unsigned int> vertices) const;

        /* -------------------- OTHER UTILITIES -----------------------*/
        void Swap();

    private:
        VertexMap<std::set<unsigned int>> _equal_constraints;
        VertexMap<std::set<unsigned int>> _diff_constraints;

        void _AtomicAddEqualConstraint(unsigned int vertex_a, unsigned int vertex_b);
        void _AtomicAddDiffConstraint(unsigned int vertex_a, unsigned int vertex_b);
};

#endif // CONSTRAINTS_HPP


/*
    Possible future implementation:
        A matrix of bytes. -> 1 matrix to represent all
            In the upper triangular matrix, (i,j) = 1, where i<j, means that nodes i and j must have same color
            In the lower triangular matrix, (i,j) = 1, where i>j, means that nodes i and j must have different color
        Potentially could be a matrix of BITES (I mean, a matrix of byte but each byte contains data about 8 vertices)

        Advantages:
            - compactness
            - faster to be copied
        Downsides:
            - access to all the constrains of a certain vertex is more inefficient
                - with map&set: O(s), where s is the number of constraints, and usually s << n° vertices
                - with matrix:  O(n° vertices), I have to scroll the entire row
                    - also, I would always read half a row and half a column. Reading a column is not cache friendly

    So no, the current implementation is fine
*/
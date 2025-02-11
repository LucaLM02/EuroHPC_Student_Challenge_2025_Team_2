#ifndef COMMON_HPP
#define COMMON_HPP

#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <memory>
#include "graph.hpp"

template <class Value>
using VertexMap = std::map<unsigned int, Value, std::less<unsigned int>>;

using VertexSet = std::vector<unsigned int>;  // comparator needs to be added

template <class Value>
using VertexPairMap = std::map<std::pair<unsigned int, unsigned int>,
			       Value>;	// comparator needs to be added

template <class Value>
using VertexPairSet = std::set<std::pair<unsigned int, unsigned int>,
			       Value>;	// comparator needs to be added

template <class Value>
using EdgeMap = VertexPairMap<Value>;  // comparator needs to be added

using Edges = std::vector<std::vector<bool>>;

void GetNeighbours(const Edges& edges, const unsigned int vertex_index,
		   VertexSet& neighbours);

/**
 * @brief Structure to manage branches and priorities in parallel
 * branch-and-bound solver.
 */
using GraphPtr = std::unique_ptr<Graph>;
struct Branch {
	GraphPtr g;
	int lb;
	unsigned short ub;
	int depth;
	bool operator<(const Branch& other) const {
		return depth < other.depth;
	}

	Branch() = default;
	Branch(GraphPtr graph, int lower, unsigned short upper, int dp) 
        : g(std::move(graph)), lb(lower), ub(upper), depth(dp) {}
};

#endif	// COMMON_HPP

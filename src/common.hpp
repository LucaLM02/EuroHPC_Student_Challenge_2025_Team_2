#ifndef COMMON_HPP
#define COMMON_HPP

#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>
#include <cstring>

#include "graph.hpp"
#include "csr_graph.hpp"

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
	Branch& operator=(const Branch& other) {
		g = std::move(other.g->Clone());
		lb = other.lb;
		ub = other.ub;
		depth = other.depth;
		return *this;
	}
	Branch() = default;
	Branch(const Branch& b)
	    : g(std::move(b.g->Clone())), lb(b.lb), ub(b.ub), depth(b.depth) {};
	Branch(GraphPtr graph, int lower, unsigned short upper, int dp)
	    : g(std::move(graph)), lb(lower), ub(upper), depth(dp) {}

	// Method to serialize branch
	std::vector<char> serialize() const {
		std::vector<char> buffer;
		std::string graphData = g->Serialize();
		size_t graphSize = graphData.size();
	
		buffer.resize(sizeof(lb) + sizeof(ub) + sizeof(depth) + sizeof(graphSize) + graphSize);
	
		char* ptr = buffer.data();
		std::memcpy(ptr, &lb, sizeof(lb));
		ptr += sizeof(lb);
		std::memcpy(ptr, &ub, sizeof(ub));
		ptr += sizeof(ub);
		std::memcpy(ptr, &depth, sizeof(depth));
		ptr += sizeof(depth);
		std::memcpy(ptr, &graphSize, sizeof(graphSize));
		ptr += sizeof(graphSize);
		std::memcpy(ptr, graphData.data(), graphSize);
	
		return buffer;
	}

	// Method to deserialize branch
	static Branch deserialize(const std::vector<char>& buffer) {
		Branch b;
		const char* ptr = buffer.data();
	
		std::memcpy(&b.lb, ptr, sizeof(b.lb));
		ptr += sizeof(b.lb);
		std::memcpy(&b.ub, ptr, sizeof(b.ub));
		ptr += sizeof(b.ub);
		std::memcpy(&b.depth, ptr, sizeof(b.depth));
		ptr += sizeof(b.depth);
	
		size_t graphSize;
		std::memcpy(&graphSize, ptr, sizeof(graphSize));
		ptr += sizeof(graphSize);
	
		std::string graphData(ptr, graphSize);
		// TODO: Add support to dynamically unserialize different graph types
		b.g = std::make_unique<CSRGraph>();
		b.g->Deserialize(graphData);
		
		return b;
	}
};

#endif	// COMMON_HPP

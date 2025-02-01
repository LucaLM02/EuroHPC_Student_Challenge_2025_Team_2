#ifndef BRANCH_N_BOUND_SEQ_HPP
#define BRANCH_N_BOUND_SEQ_HPP

#include <algorithm>
#include <queue>
#include <utility>
#include <vector>

#include "bound_strategies.hpp"
#include "branching_strategy.hpp"
#include "graph.hpp"

class BranchNBoundSeq {
       private:
	BranchingStrategy& _branching_strat;
	CliqueStrategy& _clique_strat;
	ColorStrategy& _color_strat;

	// Static private helper method to check timeout
	static bool CheckTimeout(
	    const std::chrono::steady_clock::time_point& start_time,
	    int timeout_seconds);

       public:
	BranchNBoundSeq(BranchingStrategy& branching_strat,
			LBStrategy& clique_strat, UBStrategy& color_strat)
	    : _branching_strat(branching_strat),
	      _clique_strat(clique_strat),
	      _color_strat(color_strat) {}

	int Solve(Graph& g);
};

#endif	// BRANCH_N_BOUND_SEQ_HPP

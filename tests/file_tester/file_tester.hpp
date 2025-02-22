#ifndef FILE_TESTER_HPP
#define FILE_TESTER_HPP

#include <iostream>
#include <filesystem>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <string>
#include <mpi.h>

#include "branching_strategy.hpp"
#include "clique_strategy.hpp"
#include "fastwclq.hpp"
#include "color.hpp"
#include "dsatur_color.hpp"
#include "csr_graph.hpp"
#include "dimacs.hpp"
#include "test_common.hpp"
#include "branch_n_bound_seq.hpp"
#include "branch_n_bound_par.hpp"

class FileTester {
private:
    std::string directory;
    std::vector<std::string> fileList;
    std::unordered_map<std::string, int> expectedResults;

public:
    void init_results();
    explicit FileTester(const std::string& dir);
    void loadFiles();
    void runTests();
};

#endif // FILE_TESTER_HPP

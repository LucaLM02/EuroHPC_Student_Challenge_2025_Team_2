#include "file_tester.hpp"

namespace fs = std::filesystem;

void FileTester::init_results() {
    expectedResults.insert({"anna.col", 11});
    expectedResults.insert({"david.col", 11});
    expectedResults.insert({"fpsol2.i.1.col", 65});
    expectedResults.insert({"fpsol2.i.2.col", 30});
    expectedResults.insert({"fpsol2.i.3.col", 30});
    expectedResults.insert({"games120.col", 9});
    expectedResults.insert({"homer.col", 13});
    expectedResults.insert({"huck.col", 11});
    expectedResults.insert({"inithx.i.1.col", 54});
    expectedResults.insert({"inithx.i.2.col", 31});
    expectedResults.insert({"inithx.i.3.col", 31});
    expectedResults.insert({"jean.col", 10});
    expectedResults.insert({"le450_5a.col", 5});
    expectedResults.insert({"le450_5b.col", 5});
    expectedResults.insert({"le450_5c.col", 5});
    expectedResults.insert({"le450_5d.col", 5});
    expectedResults.insert({"le450_15a.col", 15});
    expectedResults.insert({"le450_15b.col", 15});
    expectedResults.insert({"le450_15c.col", 15});
    expectedResults.insert({"le450_15d.col", 15});
    expectedResults.insert({"le450_25a.col", 25});
    expectedResults.insert({"le450_25b.col", 25});
    expectedResults.insert({"le450_25c.col", 25});
    expectedResults.insert({"le450_25d.col", 25});
    expectedResults.insert({"miles250.col", 8});
    expectedResults.insert({"miles500.col", 20});
    expectedResults.insert({"miles750.col", 31});
    expectedResults.insert({"miles1000.col", 42});
    expectedResults.insert({"miles1500.col", 73});
    expectedResults.insert({"mulsol.i.1.col", 49});
    expectedResults.insert({"mulsol.i.2.col", 31});
    expectedResults.insert({"mulsol.i.3.col", 31});
    expectedResults.insert({"mulsol.i.4.col", 31});
    expectedResults.insert({"mulsol.i.5.col", 31});
    expectedResults.insert({"myciel3.col", 4});
    expectedResults.insert({"myciel4.col", 5});
    expectedResults.insert({"myciel5.col", 6});
    expectedResults.insert({"myciel6.col", 7});
    expectedResults.insert({"myciel7.col", 8});
    expectedResults.insert({"queen5_5.col", 5});
    expectedResults.insert({"queen6_6.col", 7});
    expectedResults.insert({"queen7_7.col", 7});
    expectedResults.insert({"queen8_8.col", 9});
    expectedResults.insert({"queen8_12.col", 12});
    expectedResults.insert({"queen9_9.col", 10});
    expectedResults.insert({"queen11_11.col", 11});
    expectedResults.insert({"queen13_13.col", 13});
    expectedResults.insert({"zeroin.i.1.col", 49});
    expectedResults.insert({"zeroin.i.2.col", 30});
    expectedResults.insert({"zeroin.i.3.col", 30});
}

FileTester::FileTester(const std::string& dir) : directory(dir) {
    init_results();
    loadFiles();
}

void FileTester::loadFiles() {
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".col") {
            fileList.push_back(entry.path().string());
        }
    }
}

void FileTester::runTests() {
    NeighboursBranchingStrategy branching_strategy;
    FastCliqueStrategy clique_strategy;
    DSaturColorStrategy color_strategy;
    int chromatic_number;
    int expected;
    CSRGraph* graph;
    for (const auto& file : fileList) {
        std::cout << "Testing file: " << file << std::endl;

        graph = CSRGraph::LoadFromDimacs(file);

        BranchNBoundSeq solver(branching_strategy, clique_strategy, color_strategy, "log.txt");

        chromatic_number = solver.Solve(*graph, 10, 100000);

        std::string fileName = std::filesystem::path(file).filename().string();
        expected = expectedResults[fileName];

        if(chromatic_number != expected) std::cout << "Test failed: expected " << expected << " but got " << chromatic_number << std::endl;
        else std::cout << "Test passed" << std::endl;
    }
}

int main() {
    FileTester tester("graphs_instances");
    tester.runTests();
    return 0;
}


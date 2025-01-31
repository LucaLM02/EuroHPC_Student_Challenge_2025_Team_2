
#include "graph_loader.hpp"
#include <utility>
#include <vector>
#include <string>
#include <algorithm>


void GraphLoader::calculateGraphStats(int& maxDegree, int& minDegree, std::vector<float>& degreeHistogram) {
    auto degCopy = getDegrees();
    std::sort(degCopy.begin(), degCopy.end());
    size_t n = degCopy.size();
    size_t m = degreeHistogram.size();
    int cnt = 0;
    size_t di = n-1;
    for (size_t i = 0; i < m; ++i) {
        float bound = (m-i-1)*(float)n / m;
        while (degCopy[di] > bound && di > 0) {
            --di;
            ++cnt;
        }
        degreeHistogram[i] = cnt / (float)n;
        cnt = 0;
    }
    maxDegree = degCopy.front();
    minDegree = degCopy.back();
}

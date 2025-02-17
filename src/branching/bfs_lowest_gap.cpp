#include "bfs_lowest_gap.hpp"
#include "graph.hpp"

#include <algorithm>
#include <set>

// Function to select two branching vertices using DSATUR heuristic
// Picks two vertices with the highest saturation degree (number of different colors among their neighbors).
// If multiple vertices have the same saturation, it picks the ones with the highest degree.
std::pair<int, int> selectBranchingVertices(const Graph &graph, const std::vector<int> &coloring)
{
    int num_vertices = graph.GetNumVertices();
    std::vector<int> vertices = graph.GetVertices();
    std::vector<int> local_coloring(graph.GetFullColoring().begin(), graph.GetFullColoring().end());
    std::vector<int> saturation(num_vertices, 0);
    std::set<int> used_colors;

    int best_vertex_1 = -1, best_vertex_2 = -1;
    int max_saturation_1 = -1, max_saturation_2 = -1;
    int max_degree_1 = -1, max_degree_2 = -1;

    // Iterate through all vertices
    for (int v : vertices) {
        if (local_coloring[v] != -1) continue; // Skip already colored vertices

        used_colors.clear();
        std::vector<int> neighbors;
        graph.GetNeighbours(v, neighbors);

        // Compute the number of different colors among neighbors
        for (int neighbor : neighbors) {
            if (local_coloring[neighbor] != -1) {
                used_colors.insert(local_coloring[neighbor]);
            }
        }

        saturation[v] = used_colors.size();

        // Select the top two vertices based on saturation and degree
        if (saturation[v] > max_saturation_1 ||
           (saturation[v] == max_saturation_1 && neighbors.size() > max_degree_1)) {
            best_vertex_2 = best_vertex_1;
            max_saturation_2 = max_saturation_1;
            max_degree_2 = max_degree_1;

            best_vertex_1 = v;
            max_saturation_1 = saturation[v];
            max_degree_1 = neighbors.size();
        } else if (saturation[v] > max_saturation_2 ||
                  (saturation[v] == max_saturation_2 && neighbors.size() > max_degree_2)) {
            best_vertex_2 = v;
            max_saturation_2 = saturation[v];
            max_degree_2 = neighbors.size();
        }
    }

    return {best_vertex_1, best_vertex_2};
}



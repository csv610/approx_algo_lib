/**
 * Chapter 3: Metric TSP (Christofides 1.5 Approximation)
 * 
 * Theory:
 *   For symmetric TSP instances satisfying the triangle inequality, Christofides'
 *   heuristic combines the Minimum Spanning Tree (MST) and a Minimum Weight
 *   Perfect Matching on the odd-degree vertices of the MST to construct a Eulerian
 *   multigraph, which is then shortcut into a TSP tour.
 *   This achieves a guaranteed 1.5-approximation ratio.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 3: Metric TSP (Christofides 1.5-approx) ===\n";

    WeightedGraph graph;
    graph[0][1] = 10.0; graph[1][0] = 10.0;
    graph[1][2] = 15.0; graph[2][1] = 15.0;
    graph[2][3] = 10.0; graph[3][2] = 10.0;
    graph[3][0] = 15.0; graph[0][3] = 15.0;
    graph[0][2] = 20.0; graph[2][0] = 20.0;
    graph[1][3] = 20.0; graph[3][1] = 20.0;

    std::cout << "  Input Complete Metric Graph: 4 vertices\n";

    auto [tour, cost] = tsp_christofides_1_5_approx(graph);

    std::cout << "  Shortcut TSP Tour Sequence: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        std::cout << tour[i] << (i + 1 == tour.size() ? "" : " -> ");
    }
    std::cout << "\n  Tour Cost: " << cost << "\n\n";
    return 0;
}

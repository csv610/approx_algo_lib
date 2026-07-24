/**
 * Chapter 1: Vertex Cover (Factor-2 Approximation)
 * 
 * Theory:
 *   For an unweighted graph, finding the minimum vertex cover is NP-hard.
 *   This algorithm finds a maximal matching M and takes both endpoints of
 *   each edge in M. The resulting cover is guaranteed to be at most 2 * OPT
 *   because any valid vertex cover must pick at least one endpoint of each
 *   edge in the matching M.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 1: Vertex Cover (Factor-2) ===\n";

    // Define a path graph P_5: 0-1-2-3-4
    Graph graph(5);
    graph[0] = {1}; graph[1] = {0, 2}; graph[2] = {1, 3}; graph[3] = {2, 4}; graph[4] = {3};

    std::cout << "  Input Graph: Path P_5 (5 vertices, 4 edges)\n";

    // Compute factor-2 approximation
    auto approx_cover = vertex_cover_approx_2(graph);
    // Compute exact optimal via brute force
    auto exact_cover = vertex_cover_exact_bruteforce(graph);

    std::cout << "  Approximate Cover: [";
    bool first = true;
    for (int v : approx_cover) {
        std::cout << (first ? "" : ", ") << v;
        first = false;
    }
    std::cout << "] (size " << approx_cover.size() << ")\n";

    std::cout << "  Optimal Cover:     [";
    first = true;
    for (int v : exact_cover) {
        std::cout << (first ? "" : ", ") << v;
        first = false;
    }
    std::cout << "] (size " << exact_cover.size() << ")\n";

    std::cout << "  Performance Ratio: " << static_cast<double>(approx_cover.size()) / exact_cover.size() 
              << " (theoretical bound: <= 2.0)\n\n";
    return 0;
}

/**
 * Chapter 30: Multicut in Trees (2-Approximation)
 * 
 * Theory:
 *   Given a tree and a set of demand pairs, find the minimum cost set of edges
 *   whose removal disconnects all pairs. Processes pairs bottom-up based on LCA depth,
 *   yielding a factor-2 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 30: Multicut in Trees ===\n";

    std::vector<Edge> edges = {{0, 1}, {0, 2}};
    std::vector<double> costs = {2.0, 4.0};
    std::vector<Edge> pairs = {{1, 2}};

    auto cut = multicut_in_trees(3, edges, costs, pairs);

    std::cout << "  Multicut Trees Edges Selected: " << cut.size() << "\n\n";
    return 0;
}

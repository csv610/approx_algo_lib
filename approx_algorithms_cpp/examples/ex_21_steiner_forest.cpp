/**
 * Chapter 21: Steiner Forest (2-Approximation)
 * 
 * Theory:
 *   Finds a minimum cost forest connecting terminal pairs. Grows dual variables
 *   on active components (primal-dual schema) and prunes redundant edges.
 *   Guarantees a factor-2 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 21: Steiner Forest (Primal-Dual) ===\n";

    std::vector<int> vertices = {0, 1, 2, 3};
    std::vector<std::tuple<int, int, double>> edges = {{0, 1, 2.0}, {1, 2, 3.0}, {2, 3, 2.0}};
    std::vector<Edge> pairs = {{0, 3}};

    auto forest = steiner_forest_primal_dual(vertices, edges, pairs);

    std::cout << "  Steiner Forest Edges Selected: " << forest.size() << "\n\n";
    return 0;
}

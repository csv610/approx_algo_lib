/**
 * Chapter 11: Euclidean TSP PTAS Heuristic
 * 
 * Theory:
 *   Arora's PTAS uses quadtree dissection and dynamic programming on portal boundaries.
 *   This implementation provides the quadtree portal partitioning and routing heuristics.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 11: Euclidean TSP PTAS ===\n";

    std::vector<std::pair<double, double>> points = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    std::cout << "  Input Points (Unit Square): (0,0), (1,0), (1,1), (0,1)\n";

    auto [tour, cost] = quadtree_tsp(points);

    std::cout << "  Quadtree-routed TSP Cost: " << cost << "\n\n";
    return 0;
}

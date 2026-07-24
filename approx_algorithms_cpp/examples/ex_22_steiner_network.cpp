/**
 * Chapter 22: Steiner Network (Jain's Iterative Rounding)
 * 
 * Theory:
 *   Jain's algorithm iteratively solves the LP relaxation, picks any edge e with
 *   fractional value x_e >= 0.5, sets x_e = 1, and resolves the LP on the residual graph.
 *   Guarantees a factor-2 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 22: Steiner Network (Iterative Rounding) ===\n";

    std::vector<Edge> edges = {{0, 1}, {1, 2}, {2, 0}};
    std::vector<double> costs = {2.0, 3.0, 4.0};
    std::map<Edge, int> r = {{{0, 2}, 1}};

    auto network = jain_iterative_rounding(3, edges, costs, r);

    std::cout << "  Iterative Rounding Network Edges: " << network.size() << "\n\n";
    return 0;
}

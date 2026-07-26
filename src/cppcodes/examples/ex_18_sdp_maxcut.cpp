/**
 * Chapter 18: Semidefinite Programming (Max-Cut)
 * 
 * Theory:
 *   Goemans and Williamson's Max-Cut algorithm relaxes the problem to a Semidefinite Program,
 *   embeds vertices on a unit sphere, and rounds using a random hyperplane.
 *   This achieves an approximation ratio >= 0.878 (GW constant).
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 18: Semidefinite Programming (Max-Cut) ===\n";

    std::vector<Edge> edges = {{0, 1}, {1, 2}, {2, 0}};
    std::vector<double> weights = {1.0, 1.0, 1.0};

    // Optimize vector embeddings
    auto vectors = optimize_max_cut_vectors(3, edges, weights);
    auto [cut, avg, val] = goemans_williamson_max_cut(3, edges, weights, vectors);

    std::cout << "  SDP Rounded Best Cut Value: " << val << " (bound: >= 0.878 * OPT)\n\n";
    return 0;
}

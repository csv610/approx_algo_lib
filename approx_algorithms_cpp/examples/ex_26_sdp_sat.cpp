/**
 * Chapter 26: Semidefinite Programming (Max 2-SAT)
 * 
 * Theory:
 *   Goemans-Williamson algorithm for Max 2-SAT embeds variables in a 3D sphere
 *   and rounds using a random cutting hyperplane, achieving an 0.878 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 26: Semidefinite Programming (Max 2-SAT) ===\n";

    std::vector<std::pair<int, int>> clauses = {{1, 2}, {-1, 2}};
    std::vector<double> weights = {1.0, 1.0};

    auto vectors = optimize_max_2sat_vectors(2, clauses, weights);
    auto [assign, avg, max_w] = goemans_williamson_max_2sat(2, clauses, weights, vectors);

    std::cout << "  Max 2-SAT Rounded satisfied weight: " << max_w << "\n\n";
    return 0;
}

/**
 * Chapter 23: Feedback Vertex Set via Primal-Dual
 * 
 * Theory:
 *   Uses primal-dual variable constraints on cycles to grow dual variables
 *   and select acyclic feedback vertex sets.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 23: Feedback Vertex Set (PD) ===\n";

    std::vector<int> vertices = {0, 1, 2};
    std::vector<Edge> edges = {{0, 1}, {1, 2}, {2, 0}};
    std::map<int, double> weights = {{0, 10.0}, {1, 10.0}, {2, 3.0}};

    auto fvs = primal_dual_fvs(vertices, edges, weights);

    std::cout << "  PD Feedback Vertex Set Size: " << fvs.size() << "\n\n";
    return 0;
}

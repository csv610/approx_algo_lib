/**
 * Chapter 5: k-Center (2-Approximation)
 * 
 * Theory:
 *   The bottleneck k-Center problem selects k vertices as centers to minimize the
 *   maximum distance from any vertex to its closest center. Parametric pruning
 *   computes threshold graphs and finds a maximal independent set, guaranteeing
 *   a 2-approximation (which is optimal unless P=NP).
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 5: k-Center (2-approx) ===\n";

    WeightedGraph g;
    g[0][1] = 2.0; g[1][0] = 2.0;
    g[1][2] = 3.0; g[2][1] = 3.0;
    g[0][2] = 4.0; g[2][0] = 4.0;

    std::cout << "  Input Complete Metric Graph: 3 vertices, k = 2\n";

    auto [centers, radius] = kcenter_parametric_pruning(g, 2);

    std::cout << "  Selected Centers: [";
    bool first = true;
    for (int c : centers) {
        std::cout << (first ? "" : ", ") << c;
        first = false;
    }
    std::cout << "]\n";
    std::cout << "  Bottleneck Radius: " << radius << "\n\n";
    return 0;
}

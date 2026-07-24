/**
 * Chapter 4: Multiway Cut (2 - 2/k Approximation)
 * 
 * Theory:
 *   Given terminals t_1 ... t_k, the Multiway Cut problem partitions vertices
 *   such that no two terminals belong to the same component, minimizing cut edge costs.
 *   The algorithm computes a minimum isolating cut for each terminal and discards
 *   the heaviest one, achieving a (2 - 2/k) approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 4: Multiway Cut (2 - 2/k) ===\n";

    WeightedGraph g;
    g[0][1] = 2.0; g[1][0] = 2.0;
    g[1][2] = 2.0; g[2][1] = 2.0;
    g[2][0] = 3.0; g[0][2] = 3.0;
    std::set<int> terminals = {0, 2};

    std::cout << "  Graph: Triangle (3 vertices), Terminals: {0, 2}\n";

    auto [cut, cost] = multiway_cut_2_2k(g, terminals);

    std::cout << "  Selected Cut Edges Count: " << cut.size() << "\n";
    std::cout << "  Cut Cost: " << cost << "\n\n";
    return 0;
}

/**
 * Chapter 19: Multiway Cut via LP Rounding
 * 
 * Theory:
 *   CKR randomized rounding embeds vertices in a simplex and cuts components using
 *   a random radius threshold, guaranteeing a 1.3438 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 19: Multiway Cut LP Rounding (CKR) ===\n";

    std::vector<Edge> edges = {{0, 1}, {1, 2}, {2, 0}};
    std::vector<double> costs = {2.0, 2.0, 2.0};
    std::vector<int> terminals = {0, 2};

    auto [d, lp_obj] = solve_multiway_cut_lp(3, edges, costs, terminals);
    auto [cut, cost] = calinescu_karloff_rabani_rounding(3, edges, costs, terminals, d);

    std::cout << "  CKR Cut cost: " << cost << "\n\n";
    return 0;
}

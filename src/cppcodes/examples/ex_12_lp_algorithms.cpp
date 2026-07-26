/**
 * Chapter 12: LP-Rounding for Set Cover
 * 
 * Theory:
 *   Formulates the Set Cover problem as an Integer Program, relaxes it to a Linear Program,
 *   solves it using Google OR-Tools, and rounds variables x_i >= 1/f (where f is the maximum
 *   frequency of any element). This yields a factor-f approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 12: LP Rounded Set Cover ===\n";

    std::set<int> universe = {0, 1, 2};
    std::map<int, std::set<int>> sets = {{0, {0, 1}}, {1, {1, 2}}};
    std::map<int, double> costs = {{0, 2.0}, {1, 2.0}};

    auto [picked, cost] = set_cover_lp_rounding(universe, sets, costs);

    std::cout << "  Set Cover LP Rounded Cost: " << cost << "\n\n";
    return 0;
}

/**
 * Chapter 2: Set Cover (Greedy H_n Approximation)
 * 
 * Theory:
 *   Given a universe U and a family of subsets S with costs c(s), the Set Cover
 *   problem seeks to cover U with minimum total cost. The greedy algorithm picks
 *   the set that minimizes cost-per-uncovered-element.
 *   This yields an H_d approximation, where d is the size of the largest set
 *   and H_d = 1 + 1/2 + ... + 1/d is the d-th Harmonic number.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 2: Set Cover (Greedy H_n) ===\n";

    std::set<int> universe = {0, 1, 2, 3, 4};
    std::map<int, std::set<int>> sets = {
        {0, {0, 1, 2}},
        {1, {2, 3}},
        {2, {3, 4}}
    };
    std::map<int, double> costs = {{0, 2.0}, {1, 1.5}, {2, 2.0}};

    std::cout << "  Universe size: " << universe.size() << "\n";
    std::cout << "  Sets available: S0(cost 2.0), S1(cost 1.5), S2(cost 2.0)\n";

    auto [picked, cost] = greedy_set_cover(universe, sets, costs);

    std::cout << "  Greedy Selected Sets: [";
    bool first = true;
    for (int s : picked) {
        std::cout << (first ? "" : ", ") << "S" << s;
        first = false;
    }
    std::cout << "]\n";
    std::cout << "  Total Cost: " << cost << "\n\n";
    return 0;
}

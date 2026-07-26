/**
 * Chapter 9: Bin Packing (APTAS)
 * 
 * Theory:
 *   Bin packing groups items of sizes <= 1 into minimum number of unit-capacity bins.
 *   The Asymptotic PTAS (APTAS) rounds small items, groups large items, solves
 *   LP relaxations, and packs remainder greedily, guaranteeing <= (1 + eps)*OPT + 1.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 9: Bin Packing APTAS ===\n";

    std::vector<double> items = {0.2, 0.5, 0.4, 0.7, 0.1, 0.3};
    double eps = 0.4;

    std::cout << "  Items: [0.2, 0.5, 0.4, 0.7, 0.1, 0.3], Epsilon: " << eps << "\n";

    auto bins = bin_packing_aptas(items, eps, 1.0);

    std::cout << "  APTAS Bins Used: " << bins.size() << "\n\n";
    return 0;
}

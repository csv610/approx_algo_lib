/**
 * Chapter 24: Facility Location (Primal-Dual)
 * 
 * Theory:
 *   Computes primal-dual connections for opening facilities and assigning clients.
 *   Guarantees a factor-3 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 24: Facility Location (Primal-Dual) ===\n";

    std::map<int, Facility> facilities = {
        {0, {10.0, {{0, 2.0}, {1, 4.0}}}},
        {1, {12.0, {{0, 4.0}, {1, 2.0}}}}
    };
    std::vector<int> clients = {0, 1};

    auto [facs, assigns, cost] = facility_location_primal_dual(facilities, clients);

    std::cout << "  PD Facilities Opened: " << facs.size() << ", Cost: " << cost << "\n\n";
    return 0;
}

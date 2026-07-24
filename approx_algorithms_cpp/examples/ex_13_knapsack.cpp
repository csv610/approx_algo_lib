/**
 * Chapter 13: Knapsack FPTAS (Vazirani Chapter 13 formulation)
 * 
 * Theory:
 *   Re-evaluates Knapsack value scaling FPTAS. Highlights rounding dynamic programming.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 13: Knapsack FPTAS (Ch 13) ===\n";

    std::vector<int> weights = {10, 20, 30};
    std::vector<int> values = {60, 100, 120};
    auto [items, value] = knapsack_fptas(weights, values, 50, 0.1);

    std::cout << "  Capacity: 50, epsilon: 0.1\n";
    std::cout << "  FPTAS Rounded Knapsack Value: " << value << "\n\n";
    return 0;
}

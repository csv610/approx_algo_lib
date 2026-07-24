/**
 * Chapter 8: Knapsack FPTAS
 * 
 * Theory:
 *   FPTAS yields a solution within (1 - eps)*OPT in polynomial time O(n^3 / eps).
 *   It scales values by K = (eps * V_max)/n and solves exact DP on scaled values.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 8: Knapsack FPTAS ===\n";

    std::vector<int> weights = {10, 20, 30};
    std::vector<int> values = {60, 100, 120};
    int capacity = 50;
    double epsilon = 0.1;

    auto [items_dp, val_dp] = knapsack_dp(weights, values, capacity);
    auto [items, val] = knapsack_fptas(weights, values, capacity, epsilon);

    std::cout << "  Capacity: " << capacity << ", epsilon: " << epsilon << "\n";
    std::cout << "  Exact DP Value: " << val_dp << "\n";
    std::cout << "  FPTAS Value:    " << val << " (ratio: " << static_cast<double>(val)/val_dp << ")\n\n";
    return 0;
}

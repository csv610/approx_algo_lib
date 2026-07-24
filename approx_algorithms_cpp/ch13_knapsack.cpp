#include <iostream>
#include <vector>
#include "chapters.hpp"
#include "approx_algorithms.hpp"

void demo_knapsack_ch13() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 13: Knapsack FPTAS (Reflected from Chapter 8)\n";
    std::cout << "============================================================\n";

    std::vector<int> weights = {10, 20, 30, 40, 50};
    std::vector<int> values = {60, 100, 120, 200, 250};
    int capacity = 100;

    std::cout << "\nInstance: n=" << weights.size() << ", capacity=" << capacity << "\n";

    // Call exact DP
    auto [sel_exact, val_exact] = knapsack_dp(weights, values, capacity);
    std::cout << "\nExact DP Value: " << val_exact << "\n";

    // Call FPTAS for different epsilon
    for (double eps : {0.5, 0.25, 0.1, 0.05, 0.01}) {
        auto [sel, val] = knapsack_fptas(weights, values, capacity, eps);
        double ratio = val_exact > 0 ? static_cast<double>(val) / val_exact : 0.0;
        std::cout << "  FPTAS eps=" << eps << ": value=" << val 
                  << ", ratio=" << ratio << " (bound=" << 1.0 - eps << ")\n";
    }
}

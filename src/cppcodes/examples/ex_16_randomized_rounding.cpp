/**
 * Chapter 16: Randomized Rounding (Max-SAT)
 * 
 * Theory:
 *   Solves the LP relaxation for Max-SAT and interprets the variables y_i* as probabilities.
 *   Variables are randomized: x_i is set to True with probability y_i*, guaranteeing
 *   a (1 - 1/e) approx (approx 0.632) in expectation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 16: Randomized Rounding (Max-SAT) ===\n";

    std::vector<std::pair<std::set<int>, std::set<int>>> clauses = {
        {{1, 2}, {}},    // x_1 or x_2
        {{2}, {1}},      // not x_1 or x_2
        {{1}, {2}}       // x_1 or not x_2
    };
    std::vector<double> weights = {1.0, 1.0, 1.0};
    std::vector<double> y_lp = {0.6, 0.7};

    auto [assign, avg_w, max_w] = randomized_rounding_max_sat(2, clauses, weights, y_lp);

    std::cout << "  Trials Run: 500\n";
    std::cout << "  Average cut weight satisfied: " << avg_w << "\n\n";
    return 0;
}

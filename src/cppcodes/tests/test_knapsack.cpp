#include <gtest/gtest.h>
#include "../approx_algorithms.hpp"

TEST(KnapsackTest, ExactAndFPTAS) {
    using namespace aal;
    std::vector<int> weights = {10, 20, 30};
    std::vector<int> values = {60, 100, 120};
    int capacity = 50;

    auto [items_dp, val_dp] = knapsack_dp(weights, values, capacity);
    EXPECT_EQ(val_dp, 220); // Picks item 1 (value 100) and item 2 (value 120)

    double eps = 0.1;
    auto [items_fptas, val_fptas] = knapsack_fptas(weights, values, capacity, eps);

    // Verify FPTAS guarantee: val_fptas >= (1 - eps) * val_dp
    EXPECT_GE(val_fptas, (1.0 - eps) * val_dp);
}

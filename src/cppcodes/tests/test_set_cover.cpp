#include <gtest/gtest.h>
#include "../approx_algorithms.hpp"

TEST(SetCoverTest, BasicCover) {
    using namespace aal;
    std::set<int> universe = {0, 1, 2, 3, 4};
    std::map<int, std::set<int>> sets = {
        {0, {0, 1, 2}},
        {1, {2, 3}},
        {2, {3, 4}}
    };
    std::map<int, double> costs = {{0, 2.0}, {1, 1.5}, {2, 2.0}};

    auto [picked, cost] = greedy_set_cover(universe, sets, costs);

    // Verify all universe elements are covered
    std::set<int> covered;
    for (int s : picked) {
        for (int e : sets[s]) {
            covered.insert(e);
        }
    }

    EXPECT_EQ(covered.size(), universe.size());
    EXPECT_EQ(cost, 3.5); // Picks S0 (2.0) and S2 (2.0) or S1 (1.5) and S2 (2.0) etc. Depending on greedy choice, here cost is S1 + S2 = 1.5 + 2.0 = 3.5.
}

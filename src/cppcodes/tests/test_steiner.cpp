#include <gtest/gtest.h>
#include "../approx_algorithms.hpp"

TEST(SteinerTest, MetricTSP) {
    using namespace aal;
    WeightedGraph graph;
    graph[0][1] = 10.0; graph[1][0] = 10.0;
    graph[1][2] = 15.0; graph[2][1] = 15.0;
    graph[2][3] = 10.0; graph[3][2] = 10.0;
    graph[3][0] = 15.0; graph[0][3] = 15.0;
    graph[0][2] = 20.0; graph[2][0] = 20.0;
    graph[1][3] = 20.0; graph[3][1] = 20.0;

    auto [tour, cost] = tsp_christofides_1_5_approx(graph);

    // Verify all vertices are visited
    std::set<int> visited(tour.begin(), tour.end());
    EXPECT_EQ(visited.size(), 4);
    EXPECT_LE(cost, 50.0); // Simple upper bound check
}

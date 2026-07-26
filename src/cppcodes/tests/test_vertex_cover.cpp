#include <gtest/gtest.h>
#include "../approx_algorithms.hpp"

TEST(VertexCoverTest, EmptyGraph) {
    using namespace aal;
    Graph g(0);
    auto cover = vertex_cover_approx_2(g);
    EXPECT_TRUE(cover.empty());
}

TEST(VertexCoverTest, Path5approx) {
    using namespace aal;
    // Path P_5: 0-1-2-3-4
    Graph g(5);
    g[0] = {1}; g[1] = {0, 2}; g[2] = {1, 3}; g[3] = {2, 4}; g[4] = {3};

    auto approx = vertex_cover_approx_2(g);
    auto exact = vertex_cover_exact_bruteforce(g);

    // Verify valid cover properties
    auto is_valid_cover = [](const Graph& graph, const std::set<int>& cover_set) -> bool {
        for (int u = 0; u < static_cast<int>(graph.size()); ++u) {
            for (int v : graph[u]) {
                if (cover_set.count(u) == 0 && cover_set.count(v) == 0) {
                    return false;
                }
            }
        }
        return true;
    };

    EXPECT_TRUE(is_valid_cover(g, approx));
    EXPECT_TRUE(is_valid_cover(g, exact));

    // Verify factor-2 bound: |approx| <= 2 * |exact|
    EXPECT_LE(approx.size(), 2 * exact.size());
}

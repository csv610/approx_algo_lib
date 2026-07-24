#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <print>

namespace aal {

std::pair<std::set<int>, std::map<std::pair<int, int>, double>> vertex_cover_primal_dual(
    const std::vector<int>& vertices,
    const std::vector<std::pair<int, int>>& edges,
    const std::map<int, double>& weights
) {
    std::map<std::pair<int, int>, double> y;
    for (auto [u, v] : edges) {
        auto e = std::make_pair(std::min(u, v), std::max(u, v));
        y[e] = 0.0;
    }

    std::map<int, double> vertex_dual_sums;
    for (int v : vertices) {
        vertex_dual_sums[v] = 0.0;
    }

    std::set<int> cover;

    for (auto [u, v] : edges) {
        auto e = std::make_pair(std::min(u, v), std::max(u, v));
        if (!cover.contains(u) && !cover.contains(v)) {
            double slack_u = weights.at(u) - vertex_dual_sums.at(u);
            double slack_v = weights.at(v) - vertex_dual_sums.at(v);
            double raise_amount = std::min(slack_u, slack_v);

            y[e] += raise_amount;
            vertex_dual_sums[u] += raise_amount;
            vertex_dual_sums[v] += raise_amount;

            if (std::abs(vertex_dual_sums[u] - weights.at(u)) < 1e-9) {
                cover.insert(u);
            }
            if (std::abs(vertex_dual_sums[v] - weights.at(v)) < 1e-9) {
                cover.insert(v);
            }
        }
    }
    return {cover, y};
}

} // namespace aal

using namespace aal;

void demo_weighted_vertex_cover_pd() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 15: Weighted Vertex Cover via Primal-Dual\n");
    std::print("{:=^60}\n", "");
    
    std::vector<int> vertices1 = {0, 1, 2, 3};
    std::vector<std::pair<int, int>> edges1 = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    std::map<int, double> weights1 = {{0, 3.0}, {1, 2.0}, {2, 4.0}, {3, 1.5}};
    
    std::print("\n1. Input Graph: 4-cycle C_4\n");
    std::print("  Weights: {{0: 3.0, 1: 2.0, 2: 4.0, 3: 1.5}}\n");
    
    auto [cover, y] = vertex_cover_primal_dual(vertices1, edges1, weights1);
    
    double primal_cost = 0;
    for (int v : cover) primal_cost += weights1.at(v);
    
    double dual_cost = 0;
    for (auto& [k, val] : y) dual_cost += val;
    
    std::print("  Primal Cover:    {{");
    bool first = true;
    for (int v : cover) {
        if (!first) std::print(", ");
        std::print("{}", v);
        first = false;
    }
    std::print("}} (cost: {:.2f})\n", primal_cost);
    
    std::print("  Dual Variables:  {{");
    first = true;
    for (auto& [k, val] : y) {
        if (!first) std::print(", ");
        std::print("({}, {}): {}", k.first, k.second, val);
        first = false;
    }
    std::print("}} (dual obj: {:.2f})\n", dual_cost);
    std::print("  Primal/Dual Ratio: {:.4f} (theoretical bound <= 2.00)\n", primal_cost/dual_cost);
    
    std::vector<int> vertices2 = {0, 1, 2, 3};
    std::vector<std::pair<int, int>> edges2 = {{0, 1}, {0, 2}, {0, 3}};
    std::map<int, double> weights2 = {{0, 10.0}, {1, 4.0}, {2, 4.0}, {3, 4.0}};
    
    std::print("\n2. Star Graph (center weight 10, leaves weight 4)\n");
    auto [cover2, y2] = vertex_cover_primal_dual(vertices2, edges2, weights2);
    
    double primal_cost2 = 0;
    for (int v : cover2) primal_cost2 += weights2.at(v);
    
    double dual_cost2 = 0;
    for (auto& [k, val] : y2) dual_cost2 += val;
    
    std::print("  Primal Cover:    {{");
    first = true;
    for (int v : cover2) {
        if (!first) std::print(", ");
        std::print("{}", v);
        first = false;
    }
    std::print("}} (cost: {:.2f})\n", primal_cost2);
    
    std::print("  Dual Variables:  {{");
    first = true;
    for (auto& [k, val] : y2) {
        if (!first) std::print(", ");
        std::print("({}, {}): {}", k.first, k.second, val);
        first = false;
    }
    std::print("}} (dual obj: {:.2f})\n", dual_cost2);
    std::print("  Primal/Dual Ratio: {:.4f}\n", primal_cost2/dual_cost2);
}

#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <algorithm>
#include <print>

namespace aal {

std::pair<std::map<int, int>, std::map<int, int>> get_tree_properties(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    int root = 0
) {
    std::map<int, std::vector<int>> adj;
    for (int i = 0; i < n; ++i) adj[i] = {};
    for (auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    
    std::map<int, int> parent;
    std::map<int, int> depth;
    parent[root] = -1;
    depth[root] = 0;
    
    std::vector<int> queue = {root};
    int head = 0;
    while (head < queue.size()) {
        int curr = queue[head++];
        for (int neighbor : adj[curr]) {
            if (neighbor != parent[curr]) {
                parent[neighbor] = curr;
                depth[neighbor] = depth[curr] + 1;
                queue.push_back(neighbor);
            }
        }
    }
    
    return {parent, depth};
}

int get_lca(int u, int v, const std::map<int, int>& parent, const std::map<int, int>& depth) {
    int u_curr = u;
    int v_curr = v;
    while (depth.at(u_curr) > depth.at(v_curr)) {
        u_curr = parent.at(u_curr);
    }
    while (depth.at(v_curr) > depth.at(u_curr)) {
        v_curr = parent.at(v_curr);
    }
    while (u_curr != v_curr) {
        u_curr = parent.at(u_curr);
        v_curr = parent.at(v_curr);
    }
    return u_curr;
}

std::vector<std::pair<int, int>> get_path_edges(int u, int v, const std::map<int, int>& parent) {
    std::vector<int> path_u;
    int curr = u;
    while (curr != -1) {
        path_u.push_back(curr);
        curr = parent.at(curr);
    }
    
    std::vector<int> path_v;
    curr = v;
    while (curr != -1) {
        path_v.push_back(curr);
        curr = parent.at(curr);
    }
    
    int lca_node = -1;
    std::set<int> set_v(path_v.begin(), path_v.end());
    for (int node : path_u) {
        if (set_v.contains(node)) {
            lca_node = node;
            break;
        }
    }
    
    std::vector<std::pair<int, int>> edges_on_path;
    curr = u;
    while (curr != lca_node) {
        int p = parent.at(curr);
        int a = std::min(curr, p);
        int b = std::max(curr, p);
        edges_on_path.push_back({a, b});
        curr = p;
    }
    
    curr = v;
    while (curr != lca_node) {
        int p = parent.at(curr);
        int a = std::min(curr, p);
        int b = std::max(curr, p);
        edges_on_path.push_back({a, b});
        curr = p;
    }
    
    return edges_on_path;
}

std::vector<std::pair<int, int>> multicut_in_trees(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& costs,
    const std::vector<std::pair<int, int>>& pairs
) {
    auto [parent, depth] = get_tree_properties(n, edges);
    
    std::vector<std::tuple<int, int, int>> pair_lcas;
    for (size_t i = 0; i < pairs.size(); ++i) {
        int u = pairs[i].first;
        int v = pairs[i].second;
        int lca = get_lca(u, v, parent, depth);
        pair_lcas.push_back({depth[lca], static_cast<int>(i), lca});
    }
    
    std::sort(pair_lcas.begin(), pair_lcas.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) > std::get<0>(b);
    });
    
    std::map<std::pair<int, int>, int> edge_to_idx;
    for (size_t idx = 0; idx < edges.size(); ++idx) {
        int a = std::min(edges[idx].first, edges[idx].second);
        int b = std::max(edges[idx].first, edges[idx].second);
        edge_to_idx[{a, b}] = idx;
    }
    
    std::vector<double> load(edges.size(), 0.0);
    std::vector<std::pair<int, int>> chosen_edges;
    
    for (auto& tup : pair_lcas) {
        int pair_idx = std::get<1>(tup);
        int u = pairs[pair_idx].first;
        int v = pairs[pair_idx].second;
        auto path_edges = get_path_edges(u, v, parent);
        
        bool already_cut = false;
        for (auto& e : path_edges) {
            if (std::find(chosen_edges.begin(), chosen_edges.end(), e) != chosen_edges.end()) {
                already_cut = true;
                break;
            }
        }
        
        if (already_cut) continue;
        
        double min_slack = std::numeric_limits<double>::infinity();
        std::pair<int, int> best_edge = {-1, -1};
        bool found = false;
        
        for (auto& e : path_edges) {
            int e_idx = edge_to_idx[e];
            double slack = costs[e_idx] - load[e_idx];
            if (slack < min_slack - 1e-9) {
                min_slack = slack;
                best_edge = e;
                found = true;
            } else if (std::abs(slack - min_slack) < 1e-9) {
                int d1 = std::min(depth[e.first], depth[e.second]);
                int d2 = std::min(depth[best_edge.first], depth[best_edge.second]);
                if (d1 < d2) {
                    best_edge = e;
                }
            }
        }
        
        for (auto& e : path_edges) {
            int e_idx = edge_to_idx[e];
            load[e_idx] += min_slack;
        }
        
        if (found) {
            chosen_edges.push_back(best_edge);
        }
    }
    
    std::vector<std::pair<int, int>> pruned = chosen_edges;
    for (int i = chosen_edges.size() - 1; i >= 0; --i) {
        auto e = chosen_edges[i];
        auto it = std::find(pruned.begin(), pruned.end(), e);
        if (it != pruned.end()) {
            pruned.erase(it);
        }
        
        bool all_cut = true;
        for (auto& p : pairs) {
            auto path_edges = get_path_edges(p.first, p.second, parent);
            bool cut_found = false;
            for (auto& pe : path_edges) {
                if (std::find(pruned.begin(), pruned.end(), pe) != pruned.end()) {
                    cut_found = true;
                    break;
                }
            }
            if (!cut_found) {
                all_cut = false;
                break;
            }
        }
        
        if (!all_cut) {
            pruned.push_back(e);
        }
    }
    
    return pruned;
}

} // namespace aal

using namespace aal;

void demo_tree_multicut() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 30: Multicut in Trees\n");
    std::print("{:=^60}\n", "");
    
    int n = 7;
    std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}};
    std::vector<double> costs = {2.0, 4.0, 1.0, 3.0, 2.0, 2.0};
    
    std::vector<std::pair<int, int>> pairs = {{3, 4}, {5, 6}, {3, 5}};
    
    std::print("\n1. Tree Instance and Demand Pairs:\n");
    std::print("  Edges & Costs:\n");
    for (size_t i = 0; i < edges.size(); ++i) {
        std::print("    Edge {}: ({}, {}) cost={:.1f}\n", i, edges[i].first, edges[i].second, costs[i]);
    }
    std::print("  Pairs to disconnect:\n");
    for (size_t i = 0; i < pairs.size(); ++i) {
        std::print("    Pair {}: ({}, {})\n", i, pairs[i].first, pairs[i].second);
    }
    
    auto chosen = multicut_in_trees(n, edges, costs, pairs);
    std::print("\n2. Multicut Results:\n");
    std::print("  Selected Cut Edges: [");
    double total_cost = 0.0;
    for (size_t i = 0; i < chosen.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("({}, {})", chosen[i].first, chosen[i].second);
        
        for (size_t j = 0; j < edges.size(); ++j) {
            if ((edges[j].first == chosen[i].first && edges[j].second == chosen[i].second) ||
                (edges[j].first == chosen[i].second && edges[j].second == chosen[i].first)) {
                total_cost += costs[j];
                break;
            }
        }
    }
    std::print("]\n");
    std::print("  Total Cut Cost:     {:.2f}\n", total_cost);
}

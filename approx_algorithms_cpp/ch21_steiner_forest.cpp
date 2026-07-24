#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <print>

namespace aal {

std::vector<std::set<int>> find_components(const std::vector<int>& vertices, const std::vector<std::pair<int, int>>& forest_edges) {
    std::map<int, std::vector<int>> adj;
    for (int v : vertices) adj[v] = {};
    for (auto& edge : forest_edges) {
        adj[edge.first].push_back(edge.second);
        adj[edge.second].push_back(edge.first);
    }
    
    std::set<int> visited;
    std::vector<std::set<int>> components;
    
    for (int v : vertices) {
        if (!visited.contains(v)) {
            std::set<int> comp;
            std::vector<int> stack = {v};
            while (!stack.empty()) {
                int curr = stack.back();
                stack.pop_back();
                if (!comp.contains(curr)) {
                    comp.insert(curr);
                    visited.insert(curr);
                    for (int nxt : adj[curr]) {
                        if (!visited.contains(nxt)) {
                            stack.push_back(nxt);
                        }
                    }
                }
            }
            components.push_back(comp);
        }
    }
    return components;
}

bool is_connected(
    const std::vector<int>& vertices,
    const std::vector<std::pair<int, int>>& forest_edges,
    int s,
    int t
) {
    std::map<int, std::vector<int>> adj;
    for (int v : vertices) adj[v] = {};
    for (auto& edge : forest_edges) {
        adj[edge.first].push_back(edge.second);
        adj[edge.second].push_back(edge.first);
    }
    
    std::set<int> visited = {s};
    std::vector<int> stack = {s};
    while (!stack.empty()) {
        int curr = stack.back();
        stack.pop_back();
        if (curr == t) return true;
        for (int nxt : adj[curr]) {
            if (!visited.contains(nxt)) {
                visited.insert(nxt);
                stack.push_back(nxt);
            }
        }
    }
    return false;
}

std::vector<std::pair<int, int>> steiner_forest_primal_dual(
    const std::vector<int>& vertices,
    const std::vector<std::tuple<int, int, double>>& edges,
    const std::vector<std::pair<int, int>>& pairs
) {
    std::vector<std::pair<int, int>> chosen_edges;
    std::vector<double> L(edges.size(), 0.0);
    
    while (true) {
        auto components = find_components(vertices, chosen_edges);
        std::map<int, int> node_to_comp;
        for (size_t idx = 0; idx < components.size(); ++idx) {
            for (int node : components[idx]) {
                node_to_comp[node] = idx;
            }
        }
        
        std::set<int> active_indices;
        for (auto& p : pairs) {
            int c_s = node_to_comp[p.first];
            int c_t = node_to_comp[p.second];
            if (c_s != c_t) {
                active_indices.insert(c_s);
                active_indices.insert(c_t);
            }
        }
        
        if (active_indices.empty()) {
            break;
        }
        
        double best_delta = std::numeric_limits<double>::infinity();
        std::vector<int> best_edge_indices;
        
        std::vector<int> edge_rates;
        for (size_t i = 0; i < edges.size(); ++i) {
            int u = std::get<0>(edges[i]);
            int v = std::get<1>(edges[i]);
            double cost = std::get<2>(edges[i]);
            
            bool already_chosen = false;
            for (auto& e : chosen_edges) {
                if ((e.first == u && e.second == v) || (e.first == v && e.second == u)) {
                    already_chosen = true;
                    break;
                }
            }
            if (already_chosen) {
                edge_rates.push_back(0);
                continue;
            }
            
            int c_u = node_to_comp[u];
            int c_v = node_to_comp[v];
            
            if (c_u == c_v) {
                edge_rates.push_back(0);
                continue;
            }
            
            int rate = 0;
            if (active_indices.contains(c_u)) rate += 1;
            if (active_indices.contains(c_v)) rate += 1;
            
            edge_rates.push_back(rate);
            if (rate > 0) {
                double delta = (cost - L[i]) / rate;
                if (delta < best_delta - 1e-9) {
                    best_delta = delta;
                    best_edge_indices = {static_cast<int>(i)};
                } else if (std::abs(delta - best_delta) < 1e-9) {
                    best_edge_indices.push_back(i);
                }
            }
        }
        
        for (size_t i = 0; i < edges.size(); ++i) {
            if (edge_rates[i] > 0) {
                L[i] += best_delta * edge_rates[i];
            }
        }
        
        for (int idx : best_edge_indices) {
            chosen_edges.push_back({std::get<0>(edges[idx]), std::get<1>(edges[idx])});
        }
    }
    
    std::vector<std::pair<int, int>> pruned_edges = chosen_edges;
    for (int i = chosen_edges.size() - 1; i >= 0; --i) {
        auto edge = chosen_edges[i];
        auto it = std::find(pruned_edges.begin(), pruned_edges.end(), edge);
        if (it != pruned_edges.end()) {
            pruned_edges.erase(it);
        }
        
        bool still_connected = true;
        for (auto& p : pairs) {
            if (!is_connected(vertices, pruned_edges, p.first, p.second)) {
                still_connected = false;
                break;
            }
        }
        if (!still_connected) {
            pruned_edges.push_back(edge);
        }
    }
    
    return pruned_edges;
}

} // namespace aal

using namespace aal;

void demo_steiner_forest() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 21: Steiner Forest via Primal-Dual\n");
    std::print("{:=^60}\n", "");
    
    std::vector<int> vertices = {0, 1, 2, 3, 4};
    std::vector<std::tuple<int, int, double>> edges = {
        {0, 1, 1.0},
        {1, 2, 2.0},
        {2, 3, 1.0},
        {3, 4, 3.0},
        {0, 4, 10.0}
    };
    
    std::vector<std::pair<int, int>> pairs = {{0, 2}, {2, 4}};
    
    std::print("\nGraph Vertices: [0, 1, 2, 3, 4]\n");
    std::print("Graph Edges:    [");
    for (size_t i = 0; i < edges.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("({}, {}, {:.1f})", std::get<0>(edges[i]), std::get<1>(edges[i]), std::get<2>(edges[i]));
    }
    std::print("]\n");
    std::print("Terminal Pairs: [(0, 2), (2, 4)]\n");
    
    auto forest = steiner_forest_primal_dual(vertices, edges, pairs);
    
    double forest_cost = 0.0;
    for (auto& e : edges) {
        bool in_forest = false;
        for (auto& fe : forest) {
            if ((fe.first == std::get<0>(e) && fe.second == std::get<1>(e)) ||
                (fe.first == std::get<1>(e) && fe.second == std::get<0>(e))) {
                in_forest = true;
                break;
            }
        }
        if (in_forest) {
            forest_cost += std::get<2>(e);
        }
    }
    
    std::print("\nSelected Steiner Forest Edges: [");
    for (size_t i = 0; i < forest.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("({}, {})", forest[i].first, forest[i].second);
    }
    std::print("]\n");
    std::print("Forest Total Cost:             {:.2f}\n", forest_cost);
}

#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <print>

namespace aal {

std::vector<int> find_any_cycle(const std::set<int>& vertices, const std::vector<std::pair<int, int>>& edges) {
    std::map<int, std::vector<int>> adj;
    for (int v : vertices) adj[v] = {};
    for (auto& e : edges) {
        if (adj.contains(e.first) && adj.contains(e.second)) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
    }
    
    std::set<int> visited;
    std::map<int, int> parent;
    
    auto dfs = [&](auto& self, int node, int p) -> std::vector<int> {
        visited.insert(node);
        parent[node] = p;
        for (int neighbor : adj[node]) {
            if (neighbor != p) {
                if (visited.contains(neighbor)) {
                    std::vector<int> cycle;
                    int curr = node;
                    while (curr != neighbor) {
                        cycle.push_back(curr);
                        curr = parent[curr];
                    }
                    cycle.push_back(neighbor);
                    return cycle;
                } else {
                    auto res = self(self, neighbor, node);
                    if (!res.empty()) return res;
                }
            }
        }
        return {};
    };
    
    for (int start : vertices) {
        if (!visited.contains(start)) {
            auto res = dfs(dfs, start, -1);
            if (!res.empty()) return res;
        }
    }
    return {};
}

bool is_acyclic(const std::set<int>& vertices, const std::vector<std::pair<int, int>>& edges, const std::set<int>& removed) {
    std::set<int> remaining_v;
    for (int v : vertices) {
        if (!removed.contains(v)) {
            remaining_v.insert(v);
        }
    }
    
    std::vector<std::pair<int, int>> remaining_e;
    for (auto& e : edges) {
        if (remaining_v.contains(e.first) && remaining_v.contains(e.second)) {
            remaining_e.push_back(e);
        }
    }
    
    return find_any_cycle(remaining_v, remaining_e).empty();
}

std::vector<int> primal_dual_fvs(
    const std::vector<int>& vertices,
    const std::vector<std::pair<int, int>>& edges,
    const std::map<int, double>& weights
) {
    std::map<int, double> w = weights;
    std::vector<int> S;
    
    std::set<int> active_vertices(vertices.begin(), vertices.end());
    std::vector<std::pair<int, int>> active_edges = edges;
    
    while (true) {
        while (true) {
            std::map<int, int> deg;
            for (int v : active_vertices) deg[v] = 0;
            for (auto& e : active_edges) {
                deg[e.first] += 1;
                deg[e.second] += 1;
            }
            std::vector<int> to_remove;
            for (int v : active_vertices) {
                if (deg[v] <= 1) {
                    to_remove.push_back(v);
                }
            }
            if (to_remove.empty()) break;
            for (int v : to_remove) {
                active_vertices.erase(v);
                std::vector<std::pair<int, int>> new_edges;
                for (auto& e : active_edges) {
                    if (e.first != v && e.second != v) {
                        new_edges.push_back(e);
                    }
                }
                active_edges = new_edges;
            }
        }
        
        if (active_vertices.empty()) break;
        
        auto cycle = find_any_cycle(active_vertices, active_edges);
        if (cycle.empty()) break;
        
        std::map<int, int> deg;
        for (int v : active_vertices) deg[v] = 0;
        for (auto& e : active_edges) {
            deg[e.first] += 1;
            deg[e.second] += 1;
        }
        
        double delta = std::numeric_limits<double>::infinity();
        for (int v : cycle) {
            int rate = deg[v] - 1;
            if (rate > 0) {
                double val = w[v] / rate;
                if (val < delta) {
                    delta = val;
                }
            }
        }
        
        int tight_vertex = -1;
        for (int v : cycle) {
            int rate = deg[v] - 1;
            w[v] -= delta * rate;
            if (w[v] <= 1e-9 && tight_vertex == -1) {
                tight_vertex = v;
            }
        }
        
        if (tight_vertex == -1) tight_vertex = cycle[0];
        
        S.push_back(tight_vertex);
        active_vertices.erase(tight_vertex);
        std::vector<std::pair<int, int>> new_edges;
        for (auto& e : active_edges) {
            if (e.first != tight_vertex && e.second != tight_vertex) {
                new_edges.push_back(e);
            }
        }
        active_edges = new_edges;
    }
    
    std::vector<int> pruned = S;
    for (int i = S.size() - 1; i >= 0; --i) {
        int v = S[i];
        auto it = std::find(pruned.begin(), pruned.end(), v);
        if (it != pruned.end()) {
            pruned.erase(it);
        }
        std::set<int> removed(pruned.begin(), pruned.end());
        std::set<int> v_set(vertices.begin(), vertices.end());
        if (!is_acyclic(v_set, edges, removed)) {
            pruned.push_back(v);
        }
    }
    
    return pruned;
}

} // namespace aal

using namespace aal;

void demo_primal_dual_fvs() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 23: Feedback Vertex Set via Primal-Dual\n");
    std::print("{:=^60}\n", "");
    
    std::vector<int> vertices = {0, 1, 2, 3, 4};
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 2}};
    std::map<int, double> weights = {{0, 10.0}, {1, 10.0}, {2, 3.0}, {3, 10.0}, {4, 10.0}};
    
    std::print("\n1. Intersecting Cycles (shared center node 2):\n");
    std::print("  Vertices: [0, 1, 2, 3, 4]\n");
    std::print("  Edges:    [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]\n");
    std::print("  Weights:  {{0: 10.0, 1: 10.0, 2: 3.0, 3: 10.0, 4: 10.0}}\n");
    
    auto fvs = primal_dual_fvs(vertices, edges, weights);
    std::print("  Selected FVS: [");
    for (size_t i = 0; i < fvs.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{}", fvs[i]);
    }
    std::print("]\n");
    
    double total_cost = 0;
    for (int v : fvs) total_cost += weights[v];
    std::print("  Total Cost:   {:.2f}\n", total_cost);
    
    std::vector<int> v2;
    for (int i = 0; i < 6; ++i) v2.push_back(i);
    std::vector<std::pair<int, int>> e2;
    for (int u : {0, 1, 2}) {
        for (int v : {3, 4, 5}) {
            e2.push_back({u, v});
        }
    }
    std::map<int, double> w2;
    for (int i = 0; i < 6; ++i) w2[i] = 1.0;
    
    std::print("\n2. Bipartite Graph K_3,3 (unit weights):\n");
    auto fvs2 = primal_dual_fvs(v2, e2, w2);
    std::print("  Selected FVS: [");
    for (size_t i = 0; i < fvs2.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{}", fvs2[i]);
    }
    std::print("]\n");
    
    double total_cost2 = 0;
    for (int v : fvs2) total_cost2 += w2[v];
    std::print("  Total Cost:   {:.2f}\n", total_cost2);
}

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace aal {

using Multigraph = std::map<int, std::vector<int>>;

int connected_components_count(const Multigraph& graph) {
    std::set<int> visited;
    int count = 0;
    for (const auto& kv : graph) {
        int v = kv.first;
        if (!visited.count(v)) {
            count++;
            std::vector<int> q = {v};
            visited.insert(v);
            while (!q.empty()) {
                int curr = q.front();
                q.erase(q.begin());
                if (graph.count(curr)) {
                    for (int nxt : graph.at(curr)) {
                        if (graph.count(nxt) && !visited.count(nxt)) {
                            visited.insert(nxt);
                            q.push_back(nxt);
                        }
                    }
                }
            }
        }
    }
    return count;
}

int cyclomatic_number(const Multigraph& graph) {
    int n = graph.size();
    if (n == 0) return 0;
    int m = 0;
    for (const auto& kv : graph) m += kv.second.size();
    m /= 2;
    int cc = connected_components_count(graph);
    return m - n + cc;
}

bool is_feedback_vertex_set(const Multigraph& graph, const std::set<int>& fvs) {
    std::set<int> remaining;
    for (const auto& kv : graph) {
        if (!fvs.count(kv.first)) remaining.insert(kv.first);
    }
    std::set<int> visited;
    
    for (int start : remaining) {
        if (!visited.count(start)) {
            std::vector<std::pair<int, int>> q = {{start, -1}};
            visited.insert(start);
            while (!q.empty()) {
                auto [curr, parent] = q.front();
                q.erase(q.begin());
                
                std::vector<int> neighbors;
                if (graph.count(curr)) {
                    for (int x : graph.at(curr)) {
                        if (remaining.count(x)) neighbors.push_back(x);
                    }
                }
                
                std::set<int> unique_neighbors(neighbors.begin(), neighbors.end());
                if (neighbors.size() != unique_neighbors.size()) return false;
                
                for (int nxt : neighbors) {
                    if (nxt == parent) continue;
                    if (visited.count(nxt)) return false;
                    visited.insert(nxt);
                    q.push_back({nxt, curr});
                }
            }
        }
    }
    return true;
}

std::set<int> feedback_vertex_set_approx(Multigraph graph, std::map<int, double> weights) {
    for (const auto& kv : graph) {
        int v = kv.first;
        if (weights.count(v) && weights[v] <= 1e-9) {
            Multigraph g_sub;
            for (const auto& [x, neighbors] : graph) {
                if (x != v) {
                    std::vector<int> new_neighbors;
                    for (int y : neighbors) {
                        if (y != v) new_neighbors.push_back(y);
                    }
                    g_sub[x] = new_neighbors;
                }
            }
            std::set<int> sub_fvs = feedback_vertex_set_approx(g_sub, weights);
            sub_fvs.insert(v);
            
            std::vector<int> to_check(sub_fvs.begin(), sub_fvs.end());
            for (int x : to_check) {
                sub_fvs.erase(x);
                if (!is_feedback_vertex_set(graph, sub_fvs)) {
                    sub_fvs.insert(x);
                }
            }
            return sub_fvs;
        }
    }
    
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<int> to_remove;
        for (const auto& [v, neighbors] : graph) {
            if (neighbors.size() <= 1) to_remove.push_back(v);
        }
        if (!to_remove.empty()) {
            for (int v : to_remove) {
                for (int nxt : graph[v]) {
                    if (graph.count(nxt)) {
                        auto& nxt_neighbors = graph[nxt];
                        nxt_neighbors.erase(std::remove(nxt_neighbors.begin(), nxt_neighbors.end(), v), nxt_neighbors.end());
                    }
                }
                graph.erase(v);
            }
            changed = true;
            continue;
        }
        
        std::vector<int> deg2_vertices;
        for (const auto& [v, neighbors] : graph) {
            if (neighbors.size() == 2) deg2_vertices.push_back(v);
        }
        for (int v : deg2_vertices) {
            if (!graph.count(v)) continue;
            int u = graph[v][0];
            int w = graph[v][1];
            if (u != w) {
                graph.erase(v);
                for (int& x : graph[u]) if (x == v) x = w;
                for (int& x : graph[w]) if (x == v) x = u;
                changed = true;
                break;
            }
        }
    }

    if (graph.empty() || cyclomatic_number(graph) == 0) return {};

    for (const auto& [v, neighbors] : graph) {
        if (std::find(neighbors.begin(), neighbors.end(), v) != neighbors.end()) {
            Multigraph g_sub;
            for (const auto& [x, n] : graph) {
                if (x != v) {
                    std::vector<int> new_n;
                    for (int y : n) if (y != v) new_n.push_back(y);
                    g_sub[x] = new_n;
                }
            }
            std::set<int> sub_fvs = feedback_vertex_set_approx(g_sub, weights);
            sub_fvs.insert(v);
            return sub_fvs;
        }
    }

    for (const auto& [u, neighbors] : graph) {
        for (int v : neighbors) {
            if (std::count(neighbors.begin(), neighbors.end(), v) > 1) {
                double eps = std::min(weights[u], weights[v]);
                std::map<int, double> weights_next = weights;
                weights_next[u] -= eps;
                weights_next[v] -= eps;
                
                std::set<int> sub_fvs = feedback_vertex_set_approx(graph, weights_next);
                if (!sub_fvs.count(u) && !sub_fvs.count(v)) {
                    if (weights_next[u] <= 1e-9) sub_fvs.insert(u);
                    else sub_fvs.insert(v);
                }
                return sub_fvs;
            }
        }
    }

    int cyc_g = cyclomatic_number(graph);
    std::map<int, int> deltas;
    for (const auto& kv : graph) {
        int v = kv.first;
        Multigraph g_without_v;
        for (const auto& [x, n] : graph) {
            if (x != v) {
                std::vector<int> new_n;
                for (int y : n) if (y != v) new_n.push_back(y);
                g_without_v[x] = new_n;
            }
        }
        deltas[v] = cyc_g - cyclomatic_number(g_without_v);
        if (deltas[v] <= 0) deltas[v] = 1;
    }

    double eps = std::numeric_limits<double>::infinity();
    for (const auto& kv : graph) {
        int v = kv.first;
        eps = std::min(eps, weights[v] / deltas[v]);
    }
    
    std::map<int, double> weights_next = weights;
    for (const auto& kv : graph) {
        int v = kv.first;
        weights_next[v] -= eps * deltas[v];
    }
    
    std::set<int> sub_fvs = feedback_vertex_set_approx(graph, weights_next);
    
    for (const auto& kv : graph) {
        int v = kv.first;
        if (weights_next[v] <= 1e-9 && !sub_fvs.count(v)) {
            sub_fvs.insert(v);
        }
    }
    
    std::vector<int> to_check(sub_fvs.begin(), sub_fvs.end());
    for (int v : to_check) {
        sub_fvs.erase(v);
        if (!is_feedback_vertex_set(graph, sub_fvs)) {
            sub_fvs.insert(v);
        }
    }
    
    return sub_fvs;
}

} // namespace aal

void demo_feedback_vertex_set() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 6: Feedback Vertex Set - 2-Approximation\n";
    std::cout << "============================================================\n";
    
    Multigraph graph = {
        {0, {3, 4, 5}}, {1, {3, 4, 5}}, {2, {3, 4, 5}},
        {3, {0, 1, 2}}, {4, {0, 1, 2}}, {5, {0, 1, 2}}
    };
    std::map<int, double> weights = {{0, 1.0}, {1, 1.0}, {2, 1.0}, {3, 1.0}, {4, 1.0}, {5, 1.0}};
    
    auto fvs = feedback_vertex_set_approx(graph, weights);
    std::cout << "\n1. Bipartite Graph K_3,3\n";
    std::cout << "  FVS size: " << fvs.size() << " (Optimal is 2)\n";
}

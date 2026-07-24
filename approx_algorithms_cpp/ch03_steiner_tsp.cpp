#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <tuple>
#include <stdexcept>

namespace aal {

using Graph = std::map<int, std::map<int, double>>;
using Edge = std::pair<int, int>;
using Path = std::vector<int>;

std::vector<Edge> mst_prim(const Graph& graph, const std::vector<int>& vertices_in = {}) {
    std::vector<int> vertices = vertices_in;
    if (vertices.empty()) {
        for (const auto& kv : graph) vertices.push_back(kv.first);
    }
    if (vertices.empty()) return {};

    int start = vertices[0];
    std::set<int> visited = {start};
    std::vector<Edge> edges;

    using State = std::tuple<double, int, int>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;

    std::set<int> v_set(vertices.begin(), vertices.end());

    for (const auto& [v, w] : graph.at(start)) {
        if (v_set.count(v)) pq.push({w, start, v});
    }

    while (!pq.empty() && visited.size() < vertices.size()) {
        auto [w, u, v] = pq.top();
        pq.pop();
        if (visited.count(v)) continue;
        visited.insert(v);
        edges.push_back({u, v});
        for (const auto& [w2, wgt] : graph.at(v)) {
            if (v_set.count(w2) && !visited.count(w2)) {
                pq.push({wgt, v, w2});
            }
        }
    }
    return edges;
}

std::map<int, double> dijkstra(const Graph& graph, int source, const std::set<int>& targets = {}) {
    std::map<int, double> dist;
    dist[source] = 0.0;
    using State = std::pair<double, int>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    pq.push({0.0, source});
    std::set<int> visited;
    
    std::set<int> target_set = targets;
    if (target_set.empty()) {
        for (const auto& kv : graph) target_set.insert(kv.first);
    }

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (visited.count(u)) continue;
        visited.insert(u);
        
        int in_targets = 0;
        for (int v : visited) {
            if (target_set.count(v)) in_targets++;
        }
        if (in_targets == target_set.size()) break;

        if (graph.count(u)) {
            for (const auto& [v, w] : graph.at(u)) {
                double nd = d + w;
                if (!dist.count(v) || nd < dist[v]) {
                    dist[v] = nd;
                    pq.push({nd, v});
                }
            }
        }
    }
    return dist;
}

Graph metric_closure(const Graph& graph, const std::set<int>& terminals) {
    Graph closure;
    for (int u : terminals) closure[u] = {};
    for (int u : terminals) {
        auto dist = dijkstra(graph, u, terminals);
        for (int v : terminals) {
            if (u != v) {
                closure[u][v] = dist[v];
                closure[v][u] = dist[v];
            }
        }
    }
    return closure;
}

std::pair<std::vector<Edge>, double> steiner_tree_2approx(const Graph& graph, const std::set<int>& terminals) {
    if (terminals.size() <= 1) return {{}, 0.0};
    
    auto closure = metric_closure(graph, terminals);
    std::vector<int> term_vec(terminals.begin(), terminals.end());
    auto mst_edges = mst_prim(closure, term_vec);
    
    std::set<Edge> all_edges;
    double total_weight = 0.0;
    
    for (const auto& [u, v] : mst_edges) {
        std::map<int, double> dist;
        std::map<int, int> prev;
        dist[u] = 0.0;
        using State = std::pair<double, int>;
        std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
        pq.push({0.0, u});
        
        while (!pq.empty()) {
            auto [d, x] = pq.top();
            pq.pop();
            if (x == v) break;
            if (d > dist[x]) continue;
            if (graph.count(x)) {
                for (const auto& [y, w] : graph.at(x)) {
                    double nd = d + w;
                    if (!dist.count(y) || nd < dist[y]) {
                        dist[y] = nd;
                        prev[y] = x;
                        pq.push({nd, y});
                    }
                }
            }
        }
        
        std::vector<int> path;
        int curr = v;
        while (curr != u) {
            path.push_back(curr);
            curr = prev[curr];
        }
        path.push_back(u);
        std::reverse(path.begin(), path.end());
        
        for (size_t i = 0; i < path.size() - 1; ++i) {
            int a = std::min(path[i], path[i+1]);
            int b = std::max(path[i], path[i+1]);
            all_edges.insert({a, b});
        }
    }
    
    std::vector<Edge> tree_edges;
    std::set<int> visited;
    
    auto dfs = [&](auto& self, int u, int p) -> void {
        visited.insert(u);
        if (graph.count(u)) {
            for (const auto& [v_adj, w] : graph.at(u)) {
                int a = std::min(u, v_adj);
                int b = std::max(u, v_adj);
                if (all_edges.count({a, b}) && !visited.count(v_adj)) {
                    tree_edges.push_back({u, v_adj});
                    self(self, v_adj, u);
                }
            }
        }
    };
    
    dfs(dfs, *terminals.begin(), -1);
    
    double final_weight = 0.0;
    for (const auto& [u, v] : tree_edges) {
        final_weight += graph.at(u).at(v);
    }
    return {tree_edges, final_weight};
}

std::pair<std::vector<int>, double> tsp_2approx_mst(const Graph& graph) {
    if (graph.size() <= 1) {
        std::vector<int> res;
        for(const auto& kv: graph) res.push_back(kv.first);
        return {res, 0.0};
    }
    
    auto mst_edges = mst_prim(graph);
    std::map<int, std::vector<int>> euler_adj;
    for(const auto& kv: graph) euler_adj[kv.first] = {};
    for (const auto& [u, v] : mst_edges) {
        euler_adj[u].push_back(v);
        euler_adj[v].push_back(u);
        euler_adj[u].push_back(v);
        euler_adj[v].push_back(u);
    }
    
    std::vector<int> stack = {graph.begin()->first};
    std::vector<int> tour;
    
    while (!stack.empty()) {
        int u = stack.back();
        if (!euler_adj[u].empty()) {
            int v = euler_adj[u].back();
            euler_adj[u].pop_back();
            auto it = std::find(euler_adj[v].begin(), euler_adj[v].end(), u);
            if (it != euler_adj[v].end()) euler_adj[v].erase(it);
            stack.push_back(v);
        } else {
            tour.push_back(stack.back());
            stack.pop_back();
        }
    }
    std::reverse(tour.begin(), tour.end());
    
    std::set<int> visited;
    std::vector<int> ham_cycle;
    for (int v : tour) {
        if (!visited.count(v)) {
            visited.insert(v);
            ham_cycle.push_back(v);
        }
    }
    ham_cycle.push_back(ham_cycle[0]);
    
    double cost = 0.0;
    for (size_t i = 0; i < ham_cycle.size() - 1; ++i) {
        cost += graph.at(ham_cycle[i]).at(ham_cycle[i+1]);
    }
    return {ham_cycle, cost};
}

std::vector<Edge> greedy_matching(const Graph& graph, const std::vector<int>& vertices) {
    std::set<int> remaining(vertices.begin(), vertices.end());
    std::vector<Edge> matching;
    
    std::vector<std::tuple<double, int, int>> edges;
    for (size_t i = 0; i < vertices.size(); ++i) {
        int u = vertices[i];
        for (size_t j = i + 1; j < vertices.size(); ++j) {
            int v = vertices[j];
            if (graph.count(u) && graph.at(u).count(v)) {
                edges.push_back({graph.at(u).at(v), u, v});
            }
        }
    }
    std::sort(edges.begin(), edges.end());
    
    for (const auto& [w, u, v] : edges) {
        if (remaining.count(u) && remaining.count(v)) {
            matching.push_back({u, v});
            remaining.erase(u);
            remaining.erase(v);
        }
    }
    return matching;
}

std::vector<Edge> min_weight_perfect_matching(const Graph& graph, const std::vector<int>& vertices) {
    int n = vertices.size();
    if (n == 0) return {};
    if (n % 2 == 1) throw std::invalid_argument("Odd number of vertices for perfect matching");
    if (n > 16) return greedy_matching(graph, vertices);
    
    std::map<int, double> dp;
    dp[0] = 0.0;
    std::map<int, std::pair<int, Edge>> parent;
    
    for (int mask = 0; mask < (1 << n); ++mask) {
        if (!dp.count(mask)) continue;
        int i = 0;
        while (i < n && (mask & (1 << i))) ++i;
        if (i >= n) continue;
        
        for (int j = i + 1; j < n; ++j) {
            if (!(mask & (1 << j))) {
                int u = vertices[i];
                int v = vertices[j];
                if (graph.count(u) && graph.at(u).count(v)) {
                    double w = graph.at(u).at(v);
                    int nmask = mask | (1 << i) | (1 << j);
                    if (!dp.count(nmask) || dp[mask] + w < dp[nmask]) {
                        dp[nmask] = dp[mask] + w;
                        parent[nmask] = {mask, {u, v}};
                    }
                }
            }
        }
    }
    
    int full = (1 << n) - 1;
    if (!dp.count(full)) return greedy_matching(graph, vertices);
    
    std::vector<Edge> matching;
    int mask = full;
    while (mask != 0) {
        auto [pmask, edge] = parent[mask];
        matching.push_back(edge);
        mask = pmask;
    }
    return matching;
}

std::pair<std::vector<int>, double> tsp_christofides_1_5_approx(const Graph& graph) {
    int n = graph.size();
    if (n <= 2) return tsp_2approx_mst(graph);
    
    auto mst_edges = mst_prim(graph);
    std::map<int, int> degree;
    for (const auto& kv : graph) degree[kv.first] = 0;
    for (const auto& [u, v] : mst_edges) {
        degree[u]++;
        degree[v]++;
    }
    std::vector<int> odd_vertices;
    for (const auto& [u, d] : degree) {
        if (d % 2 == 1) odd_vertices.push_back(u);
    }
    
    auto closure = metric_closure(graph, std::set<int>(odd_vertices.begin(), odd_vertices.end()));
    auto matching = min_weight_perfect_matching(closure, odd_vertices);
    
    std::map<int, std::vector<int>> euler_adj;
    for(const auto& kv: graph) euler_adj[kv.first] = {};
    for (const auto& [u, v] : mst_edges) {
        euler_adj[u].push_back(v);
        euler_adj[v].push_back(u);
    }
    for (const auto& [u, v] : matching) {
        euler_adj[u].push_back(v);
        euler_adj[v].push_back(u);
    }
    
    std::vector<int> stack = {graph.begin()->first};
    std::vector<int> tour;
    
    while (!stack.empty()) {
        int u = stack.back();
        if (!euler_adj[u].empty()) {
            int v = euler_adj[u].back();
            euler_adj[u].pop_back();
            auto it = std::find(euler_adj[v].begin(), euler_adj[v].end(), u);
            if (it != euler_adj[v].end()) euler_adj[v].erase(it);
            stack.push_back(v);
        } else {
            tour.push_back(stack.back());
            stack.pop_back();
        }
    }
    std::reverse(tour.begin(), tour.end());
    
    std::set<int> visited;
    std::vector<int> ham_cycle;
    for (int v : tour) {
        if (!visited.count(v)) {
            visited.insert(v);
            ham_cycle.push_back(v);
        }
    }
    ham_cycle.push_back(ham_cycle[0]);
    
    double cost = 0.0;
    for (size_t i = 0; i < ham_cycle.size() - 1; ++i) {
        cost += graph.at(ham_cycle[i]).at(ham_cycle[i+1]);
    }
    return {ham_cycle, cost};
}

std::tuple<Graph, double, double> tsp_tight_example_2approx(int n = 6) {
    double eps = 0.1;
    Graph tree;
    for (int i = 0; i < n; ++i) tree[i] = {};
    for (int i = 0; i < n - 1; ++i) {
        tree[i][i+1] = 1.0;
        tree[i+1][i] = 1.0;
    }
    for (int i = 0; i < n - 2; ++i) {
        tree[i][i+2] = 1.0 + eps;
        tree[i+2][i] = 1.0 + eps;
    }
    
    Graph graph;
    for (int i = 0; i < n; ++i) graph[i] = {};
    for (int start = 0; start < n; ++start) {
        std::vector<std::pair<int, double>> queue = {{start, 0.0}};
        std::set<int> visited = {start};
        int head = 0;
        while (head < queue.size()) {
            auto [curr, d] = queue[head++];
            graph[start][curr] = d;
            for (const auto& [nxt, w] : tree[curr]) {
                if (!visited.count(nxt)) {
                    visited.insert(nxt);
                    queue.push_back({nxt, d + w});
                }
            }
        }
    }
    
    auto [cycle, approx_cost] = tsp_2approx_mst(graph);
    double opt_cost = 2.0 + (n - 2) * (1.0 + eps);
    return {graph, approx_cost, opt_cost};
}

} // namespace aal

void demo_steiner_tsp() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 3: Steiner Tree and TSP\n";
    std::cout << "============================================================\n";
    
    std::cout << "\n1. Metric Steiner Tree (2-approx)\n";
    Graph graph = {
        {0, {{1, 1}, {2, 2}}},
        {1, {{0, 1}, {2, 1}, {3, 2}}},
        {2, {{0, 2}, {1, 1}, {4, 1}}},
        {3, {{1, 2}, {5, 1}}},
        {4, {{2, 1}, {5, 1}}},
        {5, {{3, 1}, {4, 1}}}
    };
    std::set<int> terminals = {0, 3, 5};
    auto [edges, cost] = steiner_tree_2approx(graph, terminals);
    std::cout << "  Terminals: {0, 3, 5}\n";
    std::cout << "  Cost: " << cost << "\n";
    
    std::cout << "\n2. Metric TSP 2-approx (MST double-tree)\n";
    Graph tsp_graph = {
        {0, {{1, 2}, {2, 3}, {3, 1}}},
        {1, {{0, 2}, {2, 2}, {3, 4}}},
        {2, {{0, 3}, {1, 2}, {3, 2}}},
        {3, {{0, 1}, {1, 4}, {2, 2}}}
    };
    auto [cycle1, cost1] = tsp_2approx_mst(tsp_graph);
    std::cout << "  Cost: " << cost1 << "\n";
    
    std::cout << "\n3. Metric TSP Christofides (3/2-approx)\n";
    auto [cycle2, cost2] = tsp_christofides_1_5_approx(tsp_graph);
    std::cout << "  Cost: " << cost2 << "\n";
    
    std::cout << "\n4. Tight example for 2-approx TSP\n";
    for (int n : {4, 6, 10}) {
        auto [g, approx, opt] = tsp_tight_example_2approx(n);
        std::cout << "  n=" << n << ": approx=" << approx << ", opt=" << opt << ", ratio=" << approx/opt << "\n";
    }
}

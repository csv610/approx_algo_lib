#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <queue>

namespace aal {

using Graph = std::map<int, std::map<int, double>>;
using Edge = std::pair<int, int>;

class MaxFlow {
public:
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<double>> cap;

    MaxFlow(int n) : n(n), adj(n), cap(n, std::vector<double>(n, 0.0)) {}

    void add_edge(int u, int v, double c) {
        cap[u][v] += c;
        if (std::find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) adj[u].push_back(v);
        if (std::find(adj[v].begin(), adj[v].end(), u) == adj[v].end()) adj[v].push_back(u);
    }

    double bfs(int s, int t, std::vector<int>& parent) {
        std::fill(parent.begin(), parent.end(), -1);
        parent[s] = -2;
        std::queue<std::pair<int, double>> q;
        q.push({s, std::numeric_limits<double>::infinity()});

        while (!q.empty()) {
            auto [u, flow] = q.front();
            q.pop();

            for (int v : adj[u]) {
                if (parent[v] == -1 && cap[u][v] > 1e-9) {
                    parent[v] = u;
                    double nf = std::min(flow, cap[u][v]);
                    if (v == t) return nf;
                    q.push({v, nf});
                }
            }
        }
        return 0.0;
    }

    double max_flow(int s, int t) {
        double flow = 0.0;
        std::vector<int> parent(n);
        while (true) {
            double nf = bfs(s, t, parent);
            if (nf < 1e-9) break;
            flow += nf;
            int v = t;
            while (v != s) {
                int u = parent[v];
                cap[u][v] -= nf;
                cap[v][u] += nf;
                v = u;
            }
        }
        return flow;
    }

    std::pair<std::set<int>, std::set<int>> min_cut(int s, int t) {
        max_flow(s, t);
        std::set<int> visited;
        std::vector<int> stack = {s};
        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            if (visited.count(u)) continue;
            visited.insert(u);
            for (int v : adj[u]) {
                if (cap[u][v] > 1e-9 && !visited.count(v)) {
                    stack.push_back(v);
                }
            }
        }
        std::set<int> T;
        for (int i = 0; i < n; ++i) {
            if (!visited.count(i)) T.insert(i);
        }
        return {visited, T};
    }
};

std::tuple<std::set<int>, std::set<int>, double> min_s_t_cut(const Graph& graph, int s, int t) {
    std::vector<int> vertices;
    for (const auto& kv : graph) vertices.push_back(kv.first);
    std::map<int, int> idx;
    for (size_t i = 0; i < vertices.size(); ++i) idx[vertices[i]] = i;
    int n = vertices.size();

    MaxFlow mf(n);
    for (const auto& [u, edges] : graph) {
        for (const auto& [v, c] : edges) {
            if (u < v) {
                mf.add_edge(idx[u], idx[v], c);
                mf.add_edge(idx[v], idx[u], c);
            }
        }
    }

    auto [S_idx, T_idx] = mf.min_cut(idx[s], idx[t]);
    std::set<int> S, T;
    for (int i : S_idx) S.insert(vertices[i]);
    for (int i : T_idx) T.insert(vertices[i]);

    double cut_weight = 0.0;
    for (int u : S) {
        for (int v : T) {
            if (graph.count(u) && graph.at(u).count(v)) {
                cut_weight += graph.at(u).at(v);
            }
        }
    }
    return {S, T, cut_weight};
}

std::pair<std::set<Edge>, double> multiway_cut_2_2k(const Graph& graph, const std::set<int>& terminals) {
    int k = terminals.size();
    if (k <= 2) {
        auto it = terminals.begin();
        int s = *it++;
        int t = *it;
        auto [S, T, w] = min_s_t_cut(graph, s, t);
        std::set<Edge> edges;
        for (int u : S) {
            for (int v : T) {
                if (graph.count(u) && graph.at(u).count(v)) {
                    edges.insert({std::min(u, v), std::max(u, v)});
                }
            }
        }
        return {edges, w};
    }

    std::vector<int> term_list(terminals.begin(), terminals.end());
    std::vector<std::set<Edge>> isolating_cuts;
    std::vector<double> cut_weights;

    std::vector<int> vertices;
    for (const auto& kv : graph) vertices.push_back(kv.first);
    std::map<int, int> idx;
    for (size_t i = 0; i < vertices.size(); ++i) idx[vertices[i]] = i;
    int n = vertices.size();
    int super_source = n;

    for (size_t i = 0; i < term_list.size(); ++i) {
        int s = term_list[i];
        MaxFlow mf(n + 1);
        for (const auto& [u, edges] : graph) {
            for (const auto& [v, c] : edges) {
                if (u < v) {
                    mf.add_edge(idx[u], idx[v], c);
                    mf.add_edge(idx[v], idx[u], c);
                }
            }
        }
        for (size_t j = 0; j < term_list.size(); ++j) {
            if (i != j) {
                mf.add_edge(super_source, idx[term_list[j]], std::numeric_limits<double>::infinity());
            }
        }

        mf.max_flow(super_source, idx[s]);
        
        std::set<int> visited;
        std::vector<int> stack = {super_source};
        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            if (visited.count(u)) continue;
            visited.insert(u);
            for (int v : mf.adj[u]) {
                if (mf.cap[u][v] > 1e-9 && !visited.count(v)) {
                    stack.push_back(v);
                }
            }
        }

        std::set<int> S;
        for (int v : visited) {
            if (v != super_source) S.insert(vertices[v]);
        }
        std::set<int> T;
        for (int v : vertices) {
            if (!S.count(v)) T.insert(v);
        }

        std::set<Edge> cut_edges;
        double cut_weight = 0.0;
        for (int u : S) {
            for (int v : T) {
                if (graph.count(u) && graph.at(u).count(v)) {
                    cut_edges.insert({std::min(u, v), std::max(u, v)});
                    cut_weight += graph.at(u).at(v);
                }
            }
        }
        isolating_cuts.push_back(cut_edges);
        cut_weights.push_back(cut_weight);
    }

    int max_idx = 0;
    for (size_t i = 1; i < cut_weights.size(); ++i) {
        if (cut_weights[i] > cut_weights[max_idx]) max_idx = i;
    }

    std::set<Edge> result_edges;
    double result_weight = 0.0;
    for (size_t i = 0; i < term_list.size(); ++i) {
        if (i != max_idx) {
            for (const auto& e : isolating_cuts[i]) result_edges.insert(e);
            result_weight += cut_weights[i];
        }
    }
    return {result_edges, result_weight};
}

std::pair<std::vector<Edge>, std::map<Edge, double>> gomory_hu_tree(const Graph& graph) {
    std::vector<int> vertices;
    for (const auto& kv : graph) vertices.push_back(kv.first);
    int n = vertices.size();
    if (n <= 1) return {{}, {}};

    std::map<int, int> idx;
    for (size_t i = 0; i < vertices.size(); ++i) idx[vertices[i]] = i;

    std::vector<int> parent(n, 0);
    std::vector<Edge> tree_edges;
    std::map<Edge, double> tree_weights;

    for (int i = 1; i < n; ++i) {
        int s = i;
        int t = parent[i];

        MaxFlow mf(n);
        for (const auto& [u, edges] : graph) {
            for (const auto& [v, c] : edges) {
                if (u < v) {
                    mf.add_edge(idx[u], idx[v], c);
                    mf.add_edge(idx[v], idx[u], c);
                }
            }
        }

        mf.max_flow(s, t);
        
        std::set<int> visited;
        std::vector<int> stack = {s};
        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            if (visited.count(u)) continue;
            visited.insert(u);
            for (int v : mf.adj[u]) {
                if (mf.cap[u][v] > 1e-9 && !visited.count(v)) {
                    stack.push_back(v);
                }
            }
        }

        int u_val = vertices[s];
        int v_val = vertices[t];
        double weight = 0.0;
        for (int x : visited) {
            for (int y = 0; y < n; ++y) {
                if (!visited.count(y)) {
                    int vx = vertices[x];
                    int vy = vertices[y];
                    if (graph.count(vx) && graph.at(vx).count(vy)) {
                        weight += graph.at(vx).at(vy);
                    }
                }
            }
        }

        tree_edges.push_back({u_val, v_val});
        tree_weights[{std::min(u_val, v_val), std::max(u_val, v_val)}] = weight;

        for (int j = i + 1; j < n; ++j) {
            if (parent[j] == t && visited.count(j)) {
                parent[j] = i;
            }
        }
    }
    return {tree_edges, tree_weights};
}

std::pair<std::set<Edge>, double> min_k_cut_2_2k(const Graph& graph, int k) {
    if (k <= 1) return {{}, 0.0};
    if (k >= graph.size()) {
        std::set<Edge> edges;
        double weight = 0.0;
        for (const auto& [u, adj] : graph) {
            for (const auto& [v, w] : adj) {
                if (u < v) {
                    edges.insert({u, v});
                    weight += w;
                }
            }
        }
        return {edges, weight};
    }

    auto [tree_edges, tree_weights] = gomory_hu_tree(graph);
    std::vector<Edge> sorted_edges = tree_edges;
    std::sort(sorted_edges.begin(), sorted_edges.end(), [&](const Edge& a, const Edge& b) {
        return tree_weights[{std::min(a.first, a.second), std::max(a.first, a.second)}] < 
               tree_weights[{std::min(b.first, b.second), std::max(b.first, b.second)}];
    });

    std::set<Edge> result_edges;
    double result_weight = 0.0;

    for (int i = 0; i < k - 1; ++i) {
        auto [u, v] = sorted_edges[i];
        auto [S, T, w] = min_s_t_cut(graph, u, v);
        for (int x : S) {
            for (int y : T) {
                if (graph.count(x) && graph.at(x).count(y)) {
                    result_edges.insert({std::min(x, y), std::max(x, y)});
                }
            }
        }
        result_weight += w;
    }
    return {result_edges, result_weight};
}

} // namespace aal

void demo_multiway_kcut() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 4: Multiway Cut and k-Cut\n";
    std::cout << "============================================================\n";
    
    std::cout << "\n1. Multiway Cut (2 - 2/k approx)\n";
    Graph graph = {
        {0, {{1, 1.0}, {3, 1.0}}},
        {1, {{0, 1.0}, {2, 2.0}}},
        {2, {{1, 2.0}, {3, 1.0}}},
        {3, {{2, 1.0}, {0, 1.0}}}
    };
    std::set<int> terminals = {0, 1, 2, 3};
    auto [edges, weight] = multiway_cut_2_2k(graph, terminals);
    std::cout << "  Weight: " << weight << "\n";
    std::cout << "  Approx factor: " << (2.0 - 2.0/terminals.size()) << "\n";
    
    std::cout << "\n3. Minimum k-Cut (2 - 2/k approx)\n";
    Graph graph2 = {
        {0, {{1, 1}, {2, 2}, {3, 2}}},
        {1, {{0, 1}, {2, 2}, {3, 2}}},
        {2, {{0, 2}, {1, 2}, {3, 1}}},
        {3, {{0, 2}, {1, 2}, {2, 1}}}
    };
    for (int k_val : {2, 3}) {
        auto [edges2, weight2] = min_k_cut_2_2k(graph2, k_val);
        std::cout << "  k=" << k_val << ": cut weight=" << weight2 << "\n";
    }
}

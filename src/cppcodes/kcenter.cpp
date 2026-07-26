#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace aal {

using Graph = std::map<int, std::map<int, double>>;

std::pair<std::vector<double>, std::vector<Graph>> build_graphs_by_threshold(const Graph& graph) {
    std::vector<double> edges;
    for (const auto& [u, adj] : graph) {
        for (const auto& [v, w] : adj) {
            if (u < v) edges.push_back(w);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    std::vector<Graph> graphs;
    for (double t : edges) {
        Graph g;
        for (const auto& kv : graph) g[kv.first] = {};
        for (const auto& [u, adj] : graph) {
            for (const auto& [v, w] : adj) {
                if (w <= t) g[u][v] = w;
            }
        }
        graphs.push_back(g);
    }
    return {edges, graphs};
}

Graph graph_square(const Graph& graph) {
    Graph g2;
    for (const auto& kv : graph) g2[kv.first] = {};

    for (const auto& [u, adj] : graph) {
        for (const auto& [v, w] : adj) {
            g2[u][v] = w;
        }
        for (const auto& [v, w] : adj) {
            if (graph.count(v)) {
                for (const auto& [w_adj, w2] : graph.at(v)) {
                    if (w_adj != u) {
                        double current = g2[u].count(w_adj) ? g2[u][w_adj] : std::numeric_limits<double>::infinity();
                        g2[u][w_adj] = std::min(current, w + w2);
                    }
                }
            }
        }
    }
    return g2;
}

std::set<int> maximal_independent_set(const Graph& graph) {
    std::set<int> independent;
    std::set<int> remaining;
    for (const auto& kv : graph) remaining.insert(kv.first);

    while (!remaining.empty()) {
        int v = *remaining.begin();
        independent.insert(v);
        std::set<int> to_remove = {v};
        if (graph.count(v)) {
            for (const auto& kv : graph.at(v)) to_remove.insert(kv.first);
        }
        for (int r : to_remove) remaining.erase(r);
    }
    return independent;
}

std::pair<std::set<int>, double> kcenter_parametric_pruning(const Graph& graph, int k) {
    if (k >= graph.size()) {
        std::set<int> res;
        for (const auto& kv : graph) res.insert(kv.first);
        return {res, 0.0};
    }

    auto [thresholds, graphs] = build_graphs_by_threshold(graph);

    for (size_t i = 0; i < thresholds.size(); ++i) {
        double t = thresholds[i];
        const auto& g = graphs[i];
        
        Graph g2 = graph_square(g);
        std::set<int> mis = maximal_independent_set(g2);
        
        if (mis.size() <= k) {
            std::set<int> centers = mis;
            if (centers.size() < k) {
                for (const auto& kv : graph) {
                    if (!centers.count(kv.first)) {
                        centers.insert(kv.first);
                        if (centers.size() == k) break;
                    }
                }
            }
            return {centers, t};
        }
    }

    std::set<int> res;
    for (const auto& kv : graph) res.insert(kv.first);
    return {res, thresholds.empty() ? 0.0 : thresholds.back()};
}

std::pair<std::set<int>, double> weighted_kcenter_3approx(const Graph& graph, const std::map<int, double>& weights, double W) {
    double total_w = 0.0;
    for (const auto& kv : weights) total_w += kv.second;
    if (W >= total_w) {
        std::set<int> res;
        for (const auto& kv : graph) res.insert(kv.first);
        return {res, 0.0};
    }

    auto [thresholds, graphs] = build_graphs_by_threshold(graph);

    for (size_t i = 0; i < thresholds.size(); ++i) {
        double t = thresholds[i];
        const auto& g = graphs[i];
        
        Graph g2 = graph_square(g);
        std::set<int> mis = maximal_independent_set(g2);
        
        std::set<int> selected;
        double current_weight = 0.0;
        
        for (int v : mis) {
            std::set<int> nbd = {v};
            if (g.count(v)) {
                for (const auto& kv : g.at(v)) nbd.insert(kv.first);
            }
            
            std::set<int> valid_nbd;
            for (int u : nbd) {
                if (weights.count(u)) valid_nbd.insert(u);
            }
            if (valid_nbd.empty()) continue;
            
            int best = -1;
            double min_w = std::numeric_limits<double>::infinity();
            for (int u : valid_nbd) {
                if (weights.at(u) < min_w) {
                    min_w = weights.at(u);
                    best = u;
                }
            }
            selected.insert(best);
            current_weight += min_w;
        }
        
        if (current_weight <= W) {
            return {selected, t};
        }
    }

    std::set<int> res;
    for (const auto& kv : graph) res.insert(kv.first);
    return {res, thresholds.empty() ? 0.0 : thresholds.back()};
}

} // namespace aal

void demo_kcenter() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 5: k-Center Problem\n";
    std::cout << "============================================================\n";
    
    std::cout << "\n1. Unweighted k-Center (Parametric Pruning - 2-approx)\n";
    Graph graph = {
        {0, {{1, 1}, {2, 3}}},
        {1, {{0, 1}, {2, 2}, {3, 4}}},
        {2, {{0, 3}, {1, 2}, {3, 1}, {4, 3}}},
        {3, {{1, 4}, {2, 1}, {4, 2}}},
        {4, {{2, 3}, {3, 2}}}
    };
    for (int k : {1, 2, 3}) {
        auto [centers, radius] = kcenter_parametric_pruning(graph, k);
        std::cout << "  k=" << k << ": radius=" << radius << "\n";
    }
    
    std::cout << "\n3. Weighted k-Center (3-approx)\n";
    std::map<int, double> weights = {{0, 1}, {1, 1}, {2, 1}, {3, 2}, {4, 2}, {5, 3}};
    double W = 4.0;
    auto [centers2, radius2] = weighted_kcenter_3approx(graph, weights, W);
    std::cout << "  Budget W=" << W << ", radius=" << radius2 << "\n";
}

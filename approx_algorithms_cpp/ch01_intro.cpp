#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <random>
#include <map>

namespace aal {

using Edge = std::pair<int, int>;
using Graph = std::vector<std::vector<int>>;

std::vector<Edge> maximal_matching(const Graph& graph) {
    int n = graph.size();
    std::vector<bool> matched(n, false);
    std::vector<Edge> matching;
    
    for (int u = 0; u < n; ++u) {
        if (!matched[u]) {
            for (int v : graph[u]) {
                if (!matched[v]) {
                    matching.emplace_back(u, v);
                    matched[u] = true;
                    matched[v] = true;
                    break;
                }
            }
        }
    }
    return matching;
}

std::set<int> vertex_cover_approx_2(const Graph& graph) {
    auto matching = maximal_matching(graph);
    std::set<int> cover;
    for (const auto& [u, v] : matching) {
        cover.insert(u);
        cover.insert(v);
    }
    return cover;
}

std::set<int> vertex_cover_approx_2_edge_weighted(std::vector<Edge> edges, const std::map<int, double>& weights) {
    std::vector<bool> covered(edges.size(), false);
    std::set<int> cover;
    
    std::sort(edges.begin(), edges.end(), [&weights](const Edge& a, const Edge& b) {
        return (weights.at(a.first) + weights.at(a.second)) < (weights.at(b.first) + weights.at(b.second));
    });
    
    for (size_t i = 0; i < edges.size(); ++i) {
        if (!covered[i]) {
            int u = edges[i].first;
            int v = edges[i].second;
            if (weights.at(u) < weights.at(v)) {
                cover.insert(u);
            } else {
                cover.insert(v);
            }
            for (size_t j = 0; j < edges.size(); ++j) {
                int u2 = edges[j].first;
                int v2 = edges[j].second;
                if (u == u2 || u == v2 || v == u2 || v == v2) {
                    covered[j] = true;
                }
            }
        }
    }
    return cover;
}

Graph tight_example_k_n_n(int n) {
    Graph graph(2 * n);
    for (int u = 0; u < n; ++u) {
        for (int v = n; v < 2 * n; ++v) {
            graph[u].push_back(v);
            graph[v].push_back(u);
        }
    }
    return graph;
}

std::set<int> vertex_cover_exact_bruteforce(const Graph& graph) {
    int n = graph.size();
    std::set<int> best;
    for (int i = 0; i < n; ++i) best.insert(i);
    
    for (int mask = 0; mask < (1 << n); ++mask) {
        std::set<int> cover;
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {
                cover.insert(i);
            }
        }
        bool valid = true;
        for (int u = 0; u < n; ++u) {
            for (int v : graph[u]) {
                if (!cover.contains(u) && !cover.contains(v)) {
                    valid = false;
                    break;
                }
            }
            if (!valid) break;
        }
        if (valid && cover.size() < best.size()) {
            best = cover;
        }
    }
    return best;
}

} // namespace aal

void demo_vertex_cover() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 1: Vertex Cover - Factor 2 Approximation\n";
    std::cout << "============================================================\n";
    
    std::cout << "\n1. Tight Example: K_{4,4}\n";
    auto g = tight_example_k_n_n(4);
    auto cover = vertex_cover_approx_2(g);
    auto exact = vertex_cover_exact_bruteforce(g);
    std::cout << "  Graph: K_{4,4} (8 vertices, 16 edges)\n";
    std::cout << "  Approx cover size: " << cover.size() << "\n";
    std::cout << "  Optimal cover size: " << exact.size() << "\n";
    std::cout << "  Ratio: " << (double)cover.size() / exact.size() << "\n";
    std::cout << "  (Optimal picks one side: 4 vertices)\n";
    
    std::cout << "\n2. Random Graph (10 vertices, p=0.3)\n";
    int n = 10;
    double p = 0.3;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    Graph g2(n);
    int edges_count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (dis(gen) < p) {
                g2[i].push_back(j);
                g2[j].push_back(i);
                edges_count++;
            }
        }
    }
    auto cover2 = vertex_cover_approx_2(g2);
    auto exact2 = vertex_cover_exact_bruteforce(g2);
    std::cout << "  Graph: " << n << " vertices, " << edges_count << " edges\n";
    std::cout << "  Approx cover size: " << cover2.size() << "\n";
    std::cout << "  Optimal cover size: " << exact2.size() << "\n";
    std::cout << "  Ratio: " << (double)cover2.size() / exact2.size() << "\n";
    
    std::cout << "\n3. Path Graph P_5\n";
    Graph g3(5);
    for (int i = 0; i < 4; ++i) {
        g3[i].push_back(i + 1);
        g3[i + 1].push_back(i);
    }
    auto cover3 = vertex_cover_approx_2(g3);
    auto exact3 = vertex_cover_exact_bruteforce(g3);
    std::cout << "  Approx cover: [";
    for (int v : cover3) std::cout << v << " ";
    std::cout << "]\n  Optimal cover: [";
    for (int v : exact3) std::cout << v << " ";
    std::cout << "]\n  Ratio: " << (double)cover3.size() / exact3.size() << "\n";
}

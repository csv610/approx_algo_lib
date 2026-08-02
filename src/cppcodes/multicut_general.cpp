// =============================================================================
// Multicut in General Graphs (Vazirani Chapter 20)
//
// Given a graph G=(V,E) with edge costs and pairs of vertices {(s_i, t_i)},
// find a minimum-cost set of edges whose removal disconnects every pair.
//
// Algorithms implemented:
//   1. Greedy multicut: sort edges by cost, greedily remove until all pairs
//      are disconnected. This gives an O(log k) approximation where k is the
//      number of pairs.
//   2. Shortest-path based: repeatedly find the shortest path (by cost) for
//      each still-connected pair and remove all edges on that path.
//
// Input:  JSON from stdin
// Output: JSON to stdout
// =============================================================================

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include "json.hpp"

using json = nlohmann::json;

namespace aal {

using Graph = std::map<int, std::map<int, double>>;
using Edge = std::pair<int, int>;

struct MulticutResult {
    std::vector<Edge> cut_edges;
    double cost;
    int num_disconnected;
};

bool are_connected(
    int s, int t,
    const std::set<Edge>& removed,
    const Graph& graph
) {
    if (s == t) return true;

    std::set<int> visited;
    std::vector<int> stack = {s};

    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        if (u == t) return true;
        if (visited.count(u)) continue;
        visited.insert(u);

        if (graph.count(u)) {
            for (const auto& [v, w] : graph.at(u)) {
                if (!visited.count(v)) {
                    Edge e = {std::min(u, v), std::max(u, v)};
                    if (!removed.count(e)) {
                        stack.push_back(v);
                    }
                }
            }
        }
    }
    return false;
}

bool all_pairs_disconnected(
    const std::vector<std::pair<int, int>>& pairs,
    const std::set<Edge>& removed,
    const Graph& graph
) {
    for (const auto& [s, t] : pairs) {
        if (are_connected(s, t, removed, graph)) {
            return false;
        }
    }
    return true;
}

std::vector<int> find_shortest_path(
    int s, int t,
    const std::set<Edge>& removed,
    const Graph& graph
) {
    if (s == t) return {s};

    std::map<int, double> dist;
    std::map<int, int> prev;
    using PDI = std::pair<double, int>;
    std::priority_queue<PDI, std::vector<PDI>, std::greater<PDI>> pq;

    for (const auto& [u, _] : graph) {
        dist[u] = std::numeric_limits<double>::infinity();
    }
    dist[s] = 0.0;
    pq.push({0.0, s});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;
        if (u == t) break;

        if (graph.count(u)) {
            for (const auto& [v, w] : graph.at(u)) {
                Edge e = {std::min(u, v), std::max(u, v)};
                if (removed.count(e)) continue;
                if (dist[u] + w < dist[v] - 1e-9) {
                    dist[v] = dist[u] + w;
                    prev[v] = u;
                    pq.push({dist[v], v});
                }
            }
        }
    }

    if (dist[t] >= std::numeric_limits<double>::infinity()) {
        return {};
    }

    std::vector<int> path;
    for (int v = t; v != s; v = prev[v]) {
        path.push_back(v);
    }
    path.push_back(s);
    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<Edge> path_to_edges(const std::vector<int>& path) {
    std::vector<Edge> edges;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        edges.push_back({std::min(path[i], path[i + 1]), std::max(path[i], path[i + 1])});
    }
    return edges;
}

MulticutResult greedy_multicut(
    const Graph& graph,
    const std::vector<std::pair<int, int>>& pairs,
    const std::vector<std::tuple<int, int, double>>& edge_list
) {
    std::vector<std::tuple<double, int, int>> sorted;
    for (const auto& [u, v, c] : edge_list) {
        sorted.push_back({c, u, v});
    }
    std::sort(sorted.begin(), sorted.end());

    std::set<Edge> removed;
    double total_cost = 0.0;

    for (const auto& [c, u, v] : sorted) {
        Edge e = {std::min(u, v), std::max(u, v)};
        if (removed.count(e)) continue;

        removed.insert(e);

        if (!all_pairs_disconnected(pairs, removed, graph)) {
            removed.erase(e);
            continue;
        }

        total_cost += c;
    }

    std::vector<Edge> cut_vec(removed.begin(), removed.end());

    return {cut_vec, total_cost, 0};
}

MulticutResult shortest_path_multicut(
    const Graph& graph,
    const std::vector<std::pair<int, int>>& pairs
) {
    std::set<Edge> removed;
    double total_cost = 0.0;

    for (const auto& [s, t] : pairs) {
        if (!are_connected(s, t, removed, graph)) continue;

        auto path = find_shortest_path(s, t, removed, graph);
        auto path_edges = path_to_edges(path);

        for (const auto& e : path_edges) {
            if (!removed.count(e)) {
                removed.insert(e);
                if (graph.count(e.first) && graph.at(e.first).count(e.second)) {
                    total_cost += graph.at(e.first).at(e.second);
                }
            }
        }
    }

    int disc = 0;
    for (const auto& [s, t] : pairs) {
        if (!are_connected(s, t, removed, graph)) {
            disc++;
        }
    }

    std::vector<Edge> cut_vec(removed.begin(), removed.end());
    return {cut_vec, total_cost, disc};
}

int count_disconnected(
    const std::vector<std::pair<int, int>>& pairs,
    const std::set<Edge>& removed,
    const Graph& graph
) {
    int count = 0;
    for (const auto& [s, t] : pairs) {
        if (!are_connected(s, t, removed, graph)) {
            count++;
        }
    }
    return count;
}

} // namespace aal

void solve(std::istream& in, std::ostream& out) {
    using namespace aal;

    json input;
    in >> input;

    int num_vertices = input["num_vertices"];
    Graph graph;
    for (int i = 0; i < num_vertices; ++i) {
        graph[i] = {};
    }

    std::vector<std::tuple<int, int, double>> edge_list;
    for (const auto& e : input["edges"]) {
        int u = e["u"];
        int v = e["v"];
        double cost = e["cost"];
        graph[u][v] = cost;
        graph[v][u] = cost;
        edge_list.push_back({u, v, cost});
    }

    std::vector<std::pair<int, int>> pairs;
    for (const auto& p : input["pairs"]) {
        pairs.push_back({p["s"], p["t"]});
    }

    auto greedy_result = greedy_multicut(graph, pairs, edge_list);

    int greedy_disc = 0;
    std::set<Edge> greedy_set(greedy_result.cut_edges.begin(), greedy_result.cut_edges.end());
    greedy_disc = count_disconnected(pairs, greedy_set, graph);

    auto sp_result = shortest_path_multicut(graph, pairs);

    json output;
    output["algorithm"] = "multicut_general";
    output["greedy_cost"] = greedy_result.cost;
    output["greedy_edges"] = json::array();
    for (const auto& [u, v] : greedy_result.cut_edges) {
        output["greedy_edges"].push_back({{"u", u}, {"v", v}});
    }
    output["shortest_path_cost"] = sp_result.cost;
    output["shortest_path_edges"] = json::array();
    for (const auto& [u, v] : sp_result.cut_edges) {
        output["shortest_path_edges"].push_back({{"u", u}, {"v", v}});
    }
    output["num_disconnected_greedy"] = greedy_disc;
    output["num_disconnected_sp"] = sp_result.num_disconnected;
    output["total_pairs"] = static_cast<int>(pairs.size());

    out << output.dump(2) << std::endl;
}

int main() {
    solve(std::cin, std::cout);
    return 0;
}

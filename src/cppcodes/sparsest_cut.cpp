// =============================================================================
// Sparsest Cut - Vazirani Chapter 21
// Exact brute-force enumeration of all subsets for small n
// Computes cut sparsity: c(S) / d(S, S_bar) where
//   c(S) = capacity of edges crossing cut (S, V\S)
//   d(S, S_bar) = sum of demands crossing the cut
// =============================================================================

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <numeric>
#include <sstream>

#include "json.hpp"

using json = nlohmann::json;

namespace aal {

using Graph = std::map<int, std::map<int, double>>;

struct Edge {
    int u, v;
    double capacity;
};

struct Demand {
    int u, v;
    double demand;
};

double cut_capacity(const std::set<int>& S, const Graph& graph) {
    double cap = 0.0;
    for (int u : S) {
        if (!graph.count(u)) continue;
        for (const auto& [v, w] : graph.at(u)) {
            if (!S.count(v)) {
                cap += w;
            }
        }
    }
    return cap;
}

double demand_crossing(const std::set<int>& S, const std::vector<Demand>& demands) {
    double total = 0.0;
    for (const auto& d : demands) {
        bool u_in = S.count(d.u);
        bool v_in = S.count(d.v);
        if (u_in != v_in) {
            total += d.demand;
        }
    }
    return total;
}

double sparsity(const std::set<int>& S, const Graph& graph,
                const std::vector<Demand>& demands) {
    double dc = demand_crossing(S, demands);
    if (dc < 1e-12) return std::numeric_limits<double>::infinity();
    return cut_capacity(S, graph) / dc;
}

double conductance(const std::set<int>& S, const Graph& graph) {
    double vol_S = 0.0;
    for (int u : S) {
        if (!graph.count(u)) continue;
        for (const auto& [v, w] : graph.at(u)) {
            vol_S += w;
        }
    }
    if (vol_S < 1e-12) return 0.0;
    return cut_capacity(S, graph) / vol_S;
}

std::vector<int> get_vertices(const Graph& graph,
                              const std::vector<Demand>& demands) {
    std::set<int> verts;
    for (const auto& [u, adj] : graph) {
        verts.insert(u);
        for (const auto& [v, w] : adj) verts.insert(v);
    }
    for (const auto& d : demands) {
        verts.insert(d.u);
        verts.insert(d.v);
    }
    return std::vector<int>(verts.begin(), verts.end());
}

struct SparsestCutResult {
    std::set<int> best_cut;
    double sparsest_value;
    double cut_cap;
    double demand_cross;
};

SparsestCutResult sparsest_cut_exact(const Graph& graph,
                                     const std::vector<Demand>& demands) {
    std::vector<int> vertices = get_vertices(graph, demands);
    int n = vertices.size();

    SparsestCutResult result;
    result.sparsest_value = std::numeric_limits<double>::infinity();

    for (int mask = 1; mask < (1 << n) - 1; ++mask) {
        std::set<int> S;
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {
                S.insert(vertices[i]);
            }
        }

        double s = sparsity(S, graph, demands);
        if (s < result.sparsest_value) {
            result.sparsest_value = s;
            result.best_cut = S;
            result.cut_cap = cut_capacity(S, graph);
            result.demand_cross = demand_crossing(S, demands);
        }
    }
    return result;
}

void solve(std::istream& in, std::ostream& out) {
    json input;
    in >> input;

    int num_vertices = input["num_vertices"].get<int>();

    Graph graph;
    std::vector<Edge> edges_input;
    if (input.contains("edges")) {
        for (const auto& e : input["edges"]) {
            int u = e["u"].get<int>();
            int v = e["v"].get<int>();
            double c = e["capacity"].get<double>();
            graph[u][v] = c;
            graph[v][u] = c;
            edges_input.push_back({u, v, c});
        }
    }

    std::vector<Demand> demands;
    if (input.contains("demands")) {
        for (const auto& d : input["demands"]) {
            int u = d["u"].get<int>();
            int v = d["v"].get<int>();
            double dem = d["demand"].get<double>();
            demands.push_back({u, v, dem});
        }
    }

    auto result = sparsest_cut_exact(graph, demands);

    std::vector<int> cut_vec(result.best_cut.begin(), result.best_cut.end());
    std::sort(cut_vec.begin(), cut_vec.end());

    json output;
    output["algorithm"] = "sparsest_cut";
    output["best_cut"] = cut_vec;
    output["sparsest_value"] = result.sparsest_value;
    output["cut_capacity"] = result.cut_cap;
    output["demand_crossing"] = result.demand_cross;

    out << output.dump(2) << std::endl;
}

} // namespace aal

int main() {
    aal::solve(std::cin, std::cout);
    return 0;
}

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include "simplex.hpp"

namespace aal {

std::pair<std::vector<double>, double> set_cover_lp(const std::set<int>& universe, const std::map<int, std::set<int>>& sets, const std::map<int, double>& costs) {
    std::vector<int> U(universe.begin(), universe.end());
    std::vector<int> S;
    for (const auto& pair : sets) S.push_back(pair.first);
    
    int m = S.size();
    int n = U.size();
    if (m == 0) {
        return {std::vector<double>(m, 0.0), 1e18};
    }
    
    std::vector<std::vector<double>> A_dual;
    for (int s : S) {
        std::vector<double> row;
        for (int e : U) {
            row.push_back(sets.at(s).count(e) ? 1.0 : 0.0);
        }
        A_dual.push_back(row);
    }
    
    std::vector<double> b_dual;
    for (int s : S) b_dual.push_back(costs.at(s));
    
    std::vector<double> c_dual(n, 1.0);
    
    Simplex simplex(A_dual, b_dual, c_dual);
    auto [y, opt] = simplex.solve();
    
    if (y.empty()) {
        return {std::vector<double>(m, 0.0), 1e18};
    }
    
    std::vector<double> x;
    for (int i = n; i < n + m; ++i) {
        x.push_back(simplex.obj_row[i]);
    }
    return {x, opt};
}

std::pair<std::vector<int>, double> set_cover_lp_rounding(const std::set<int>& universe, const std::map<int, std::set<int>>& sets, const std::map<int, double>& costs) {
    auto [x, opt] = set_cover_lp(universe, sets, costs);
    if (x.empty()) return {{}, 0.0};
    
    std::vector<int> S;
    for (const auto& pair : sets) S.push_back(pair.first);
    
    int freq = 0;
    for (int e : universe) {
        int cnt = 0;
        for (int s : S) {
            if (sets.at(s).count(e)) cnt++;
        }
        freq = std::max(freq, cnt);
    }
    
    if (freq == 0) return {{}, 0.0};
    
    double threshold = 1.0 / freq;
    std::vector<int> picked;
    double total_cost = 0.0;
    
    for (size_t j = 0; j < S.size(); ++j) {
        if (x[j] >= threshold - 1e-9) {
            picked.push_back(S[j]);
            total_cost += costs.at(S[j]);
        }
    }
    
    return {picked, total_cost};
}

std::pair<std::vector<int>, double> set_cover_primal_dual(const std::set<int>& universe, const std::map<int, std::set<int>>& sets, const std::map<int, double>& costs) {
    std::set<int> U = universe;
    std::set<int> covered;
    std::map<int, double> y;
    for (int e : U) y[e] = 0.0;
    
    std::vector<int> picked;
    double total_cost = 0.0;
    
    std::map<int, std::vector<int>> sets_containing;
    for (int e : U) sets_containing[e] = {};
    for (const auto& pair : sets) {
        int s = pair.first;
        for (int e : pair.second) {
            if (U.count(e)) sets_containing[e].push_back(s);
        }
    }
    
    while (covered.size() < U.size()) {
        int e = -1;
        for (int u : U) {
            if (!covered.count(u)) { e = u; break; }
        }
        
        while (true) {
            double min_slack = 1e18;
            int tight_set = -1;
            
            for (int s : sets_containing[e]) {
                double slack = costs.at(s);
                for (int ee : sets.at(s)) {
                    if (U.count(ee)) slack -= y[ee];
                }
                if (slack < min_slack) {
                    min_slack = slack;
                    tight_set = s;
                }
            }
            
            if (min_slack <= 1e-9) break;
            
            y[e] += min_slack;
            
            for (int s : sets_containing[e]) {
                double slack = costs.at(s);
                for (int ee : sets.at(s)) {
                    if (U.count(ee)) slack -= y[ee];
                }
                if (std::abs(slack) < 1e-9) {
                    tight_set = s;
                    break;
                }
            }
            
            if (tight_set != -1) break;
        }
        
        int tight_set = -1;
        for (int s : sets_containing[e]) {
            double slack = costs.at(s);
            for (int ee : sets.at(s)) {
                if (U.count(ee)) slack -= y[ee];
            }
            if (std::abs(slack) < 1e-9) {
                tight_set = s;
                break;
            }
        }
        
        if (tight_set == -1) break;
        
        if (std::find(picked.begin(), picked.end(), tight_set) == picked.end()) {
            picked.push_back(tight_set);
            total_cost += costs.at(tight_set);
            for (int el : sets.at(tight_set)) covered.insert(el);
        }
    }
    
    return {picked, total_cost};
}

std::pair<std::set<int>, double> vertex_cover_lp_rounding(const std::map<int, std::map<int, double>>& graph) {
    std::vector<std::pair<int, int>> edges;
    for (const auto& u_pair : graph) {
        int u = u_pair.first;
        for (const auto& v_pair : u_pair.second) {
            int v = v_pair.first;
            if (u < v) edges.push_back({u, v});
        }
    }
    
    std::set<int> universe;
    for (size_t i = 0; i < edges.size(); ++i) universe.insert(i);
    
    std::map<std::pair<int, int>, int> edge_to_idx;
    for (size_t i = 0; i < edges.size(); ++i) edge_to_idx[edges[i]] = i;
    
    std::vector<int> vertices;
    for (const auto& pair : graph) vertices.push_back(pair.first);
    
    std::map<int, std::set<int>> sets;
    for (int u : vertices) {
        sets[u] = {};
        for (const auto& v_pair : graph.at(u)) {
            int v = v_pair.first;
            auto e = std::make_pair(std::min(u, v), std::max(u, v));
            sets[u].insert(edge_to_idx[e]);
        }
    }
    
    std::map<int, double> costs;
    for (int v : vertices) costs[v] = 1.0;
    
    auto [x, opt] = set_cover_lp(universe, sets, costs);
    
    if (opt > 1e17) {
        return {std::set<int>(vertices.begin(), vertices.end()), 1e18};
    }
    
    std::set<int> cover;
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (x[i] >= 0.5 - 1e-9) cover.insert(vertices[i]);
    }
    
    return {cover, static_cast<double>(cover.size())};
}

} // namespace aal

using namespace aal;

void demo_lp_algorithms() {
    std::cout << "============================================================\n";
    std::cout << "Chapters 12-14: LP-Duality Based Algorithms\n";
    std::cout << "============================================================\n";
    
    std::set<int> universe = {1, 2, 3, 4, 5};
    std::map<int, std::set<int>> sets = {
        {0, {1, 2, 3}},
        {1, {3, 4, 5}},
        {2, {1, 4}},
        {3, {2, 5}}
    };
    std::map<int, double> costs = {{0, 3}, {1, 3}, {2, 2}, {3, 2}};
    
    std::cout << "\n1. Set Cover via LP Rounding\n";
    auto [x, opt] = set_cover_lp(universe, sets, costs);
    std::cout << "  LP solution: [";
    for (size_t i = 0; i < x.size(); ++i) std::cout << x[i] << (i + 1 == x.size() ? "" : ", ");
    std::cout << "]\n  LP optimal: " << opt << "\n";
    
    auto [picked, cost] = set_cover_lp_rounding(universe, sets, costs);
    std::cout << "  Rounded (f=2): [";
    for (size_t i = 0; i < picked.size(); ++i) std::cout << picked[i] << (i + 1 == picked.size() ? "" : ", ");
    std::cout << "], cost=" << cost << "\n";
    
    auto [picked2, cost2] = set_cover_primal_dual(universe, sets, costs);
    std::cout << "  Primal-Dual: [";
    for (size_t i = 0; i < picked2.size(); ++i) std::cout << picked2[i] << (i + 1 == picked2.size() ? "" : ", ");
    std::cout << "], cost=" << cost2 << "\n";
    
    std::cout << "\n2. Vertex Cover via LP Rounding (2-approx)\n";
    std::map<int, std::map<int, double>> graph = {
        {0, {{1, 1.0}, {2, 1.0}}},
        {1, {{0, 1.0}, {2, 1.0}, {3, 1.0}}},
        {2, {{0, 1.0}, {1, 1.0}, {3, 1.0}}},
        {3, {{1, 1.0}, {2, 1.0}}}
    };
    
    auto [cover, vcost] = vertex_cover_lp_rounding(graph);
    std::cout << "  Graph: C4 (cycle of 4)\n";
    std::cout << "  LP-rounded cover: {";
    int c = 0;
    for (int v : cover) {
        std::cout << v << (++c == cover.size() ? "" : ", ");
    }
    std::cout << "}, size=" << vcost << "\n";
    std::cout << "  Optimal: 2\n";
}

#include "chapters.hpp"
#include "simplex.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <print>

namespace aal {

std::vector<std::set<int>> get_all_cuts(int n) {
    std::vector<std::set<int>> cuts;
    for (int i = 1; i < (1 << (n - 1)); ++i) {
        std::set<int> cut;
        for (int j = 0; j < n; ++j) {
            if ((i & (1 << j)) > 0) {
                cut.insert(j);
            }
        }
        cuts.push_back(cut);
    }
    return cuts;
}

std::vector<double> solve_sndp_lp_phase(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& costs,
    const std::map<std::pair<int, int>, int>& r,
    const std::set<int>& fixed_edges
) {
    auto cuts = get_all_cuts(n);
    int n_cuts = cuts.size();
    int n_edges = edges.size();
    
    std::vector<int> active_edges;
    for (int i = 0; i < n_edges; ++i) {
        if (!fixed_edges.contains(i)) {
            active_edges.push_back(i);
        }
    }
    int n_active = active_edges.size();
    
    std::vector<double> f;
    for (const auto& cut : cuts) {
        int val = 0;
        for (const auto& [uv, req] : r) {
            bool u_in = cut.contains(uv.first);
            bool v_in = cut.contains(uv.second);
            if (u_in != v_in) {
                val = std::max(val, req);
            }
        }
        
        int fixed_crossing = 0;
        for (int edge_idx : fixed_edges) {
            int u = edges[edge_idx].first;
            int v = edges[edge_idx].second;
            if (cut.contains(u) != cut.contains(v)) {
                fixed_crossing += 1;
            }
        }
        
        f.push_back(std::max(0, val - fixed_crossing));
    }
    
    int n_cols = n_cuts + n_active;
    std::vector<std::vector<double>> A_dual;
    std::vector<double> b_dual;
    
    for (int j = 0; j < n_active; ++j) {
        int edge_idx = active_edges[j];
        int u = edges[edge_idx].first;
        int v = edges[edge_idx].second;
        
        std::vector<double> row(n_cols, 0.0);
        for (int s_idx = 0; s_idx < n_cuts; ++s_idx) {
            if (cuts[s_idx].contains(u) != cuts[s_idx].contains(v)) {
                row[s_idx] = 1.0;
            }
        }
        row[n_cuts + j] = -1.0;
        
        A_dual.push_back(row);
        b_dual.push_back(costs[edge_idx]);
    }
    
    std::vector<double> c_dual = f;
    for (int j = 0; j < n_active; ++j) {
        c_dual.push_back(-1.0);
    }
    
    Simplex solver(A_dual, b_dual, c_dual);
    auto [dual_sol, dual_obj] = solver.solve();
    
    if (dual_sol.empty()) {
        return std::vector<double>(n_edges, 0.0);
    }
    
    std::vector<double> x(n_edges, 0.0);
    for (int idx : fixed_edges) {
        x[idx] = 1.0;
    }
    
    for (int j = 0; j < n_active; ++j) {
        int edge_idx = active_edges[j];
        double val = solver.obj_row[n_cols + j];
        x[edge_idx] = std::max(0.0, val);
    }
    
    return x;
}

std::vector<std::pair<int, int>> jain_iterative_rounding(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& costs,
    const std::map<std::pair<int, int>, int>& r
) {
    std::set<int> fixed_edges;
    
    while (true) {
        auto x = solve_sndp_lp_phase(n, edges, costs, r, fixed_edges);
        
        auto cuts = get_all_cuts(n);
        bool requirements_satisfied = true;
        
        for (const auto& cut : cuts) {
            int req = 0;
            for (const auto& [uv, val] : r) {
                if (cut.contains(uv.first) != cut.contains(uv.second)) {
                    req = std::max(req, val);
                }
            }
            if (req == 0) continue;
            
            int capacity = 0;
            for (int idx : fixed_edges) {
                if (cut.contains(edges[idx].first) != cut.contains(edges[idx].second)) {
                    capacity += 1;
                }
            }
            if (capacity < req) {
                requirements_satisfied = false;
                break;
            }
        }
        
        if (requirements_satisfied) break;
        
        double best_val = -1.0;
        int best_idx = -1;
        for (size_t idx = 0; idx < edges.size(); ++idx) {
            if (!fixed_edges.contains(idx)) {
                if (x[idx] > best_val) {
                    best_val = x[idx];
                    best_idx = idx;
                }
            }
        }
        
        if (best_idx != -1 && best_val >= 0.4999) {
            fixed_edges.insert(best_idx);
        } else {
            if (best_idx != -1) {
                fixed_edges.insert(best_idx);
            } else {
                break;
            }
        }
    }
    
    std::vector<std::pair<int, int>> result;
    for (int idx : fixed_edges) {
        result.push_back(edges[idx]);
    }
    return result;
}

} // namespace aal

using namespace aal;

void demo_steiner_network() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 22: Steiner Network (Jain's Iterative Rounding)\n");
    std::print("{:=^60}\n", "");
    
    int n = 4;
    std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}};
    std::vector<double> costs = {2.0, 3.0, 2.0, 3.0, 4.0};
    
    std::map<std::pair<int, int>, int> r = {
        {{0, 2}, 2},
        {{1, 3}, 1}
    };
    
    std::print("\n1. Input Graph and Requirements:\n");
    std::print("  Edges & Costs:\n");
    for (size_t i = 0; i < edges.size(); ++i) {
        std::print("    Edge {}: ({}, {}) cost={:.1f}\n", i, edges[i].first, edges[i].second, costs[i]);
    }
    std::print("  Connectivity Requirements:\n");
    for (const auto& [uv, req] : r) {
        std::print("    r({}, {}) = {}\n", uv.first, uv.second, req);
    }
    
    auto chosen = jain_iterative_rounding(n, edges, costs, r);
    std::print("\n2. Iterative Rounding Result:\n");
    std::print("  Selected Edges: [");
    double total_cost = 0.0;
    for (size_t i = 0; i < chosen.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("({}, {})", chosen[i].first, chosen[i].second);
        
        for (size_t j = 0; j < edges.size(); ++j) {
            if ((edges[j].first == chosen[i].first && edges[j].second == chosen[i].second) ||
                (edges[j].first == chosen[i].second && edges[j].second == chosen[i].first)) {
                total_cost += costs[j];
                break;
            }
        }
    }
    std::print("]\n");
    std::print("  Total Network Cost: {:.2f}\n", total_cost);
}

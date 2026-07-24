#include "chapters.hpp"
#include "simplex.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <cmath>
#include <algorithm>
#include <print>

namespace aal {

std::tuple<std::vector<std::vector<double>>, double> solve_multiway_cut_lp(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& costs,
    const std::vector<int>& terminals
) {
    int k = terminals.size();
    int n_edges = edges.size();
    
    int n_vars = n * k + n_edges * k;
    
    std::vector<std::vector<double>> M;
    std::vector<double> f;
    
    for (int e = 0; e < n_edges; ++e) {
        int u = edges[e].first;
        int v = edges[e].second;
        for (int i = 0; i < k; ++i) {
            std::vector<double> row(n_vars, 0.0);
            row[n * k + e * k + i] = 1.0;
            row[u * k + i] = -1.0;
            row[v * k + i] = 1.0;
            M.push_back(row);
            f.push_back(0.0);
        }
    }
    
    for (int e = 0; e < n_edges; ++e) {
        int u = edges[e].first;
        int v = edges[e].second;
        for (int i = 0; i < k; ++i) {
            std::vector<double> row(n_vars, 0.0);
            row[n * k + e * k + i] = 1.0;
            row[u * k + i] = 1.0;
            row[v * k + i] = -1.0;
            M.push_back(row);
            f.push_back(0.0);
        }
    }
    
    for (int v = 0; v < n; ++v) {
        std::vector<double> row(n_vars, 0.0);
        for (int i = 0; i < k; ++i) row[v * k + i] = 1.0;
        M.push_back(row);
        f.push_back(1.0);
    }
    
    for (int v = 0; v < n; ++v) {
        std::vector<double> row(n_vars, 0.0);
        for (int i = 0; i < k; ++i) row[v * k + i] = -1.0;
        M.push_back(row);
        f.push_back(-1.0);
    }
    
    for (int j = 0; j < k; ++j) {
        int s_j = terminals[j];
        for (int i = 0; i < k; ++i) {
            std::vector<double> row(n_vars, 0.0);
            if (i == j) {
                row[s_j * k + i] = 1.0;
                M.push_back(row);
                f.push_back(1.0);
            } else {
                row[s_j * k + i] = -1.0;
                M.push_back(row);
                f.push_back(0.0);
            }
        }
    }
    
    std::vector<double> c_primal(n_vars, 0.0);
    for (int e = 0; e < n_edges; ++e) {
        for (int i = 0; i < k; ++i) {
            c_primal[n * k + e * k + i] = 0.5 * costs[e];
        }
    }
    
    int n_dual_vars = f.size();
    std::vector<std::vector<double>> A_dual(n_vars, std::vector<double>(n_dual_vars, 0.0));
    for (int row_idx = 0; row_idx < n_dual_vars; ++row_idx) {
        for (int col_idx = 0; col_idx < n_vars; ++col_idx) {
            A_dual[col_idx][row_idx] = M[row_idx][col_idx];
        }
    }
    
    Simplex solver(A_dual, c_primal, f);
    auto [dual_sol, dual_obj] = solver.solve();
    
    if (dual_sol.empty()) {
        std::vector<std::vector<double>> d(n, std::vector<double>(k, 1.0 / k));
        return {d, 0.0};
    }
    
    std::vector<double> primal_vals;
    for (int var_idx = 0; var_idx < n_vars; ++var_idx) {
        double val = solver.obj_row[n_dual_vars + var_idx];
        primal_vals.push_back(std::max(0.0, val));
    }
    
    std::vector<std::vector<double>> d(n, std::vector<double>(k, 0.0));
    for (int v = 0; v < n; ++v) {
        for (int i = 0; i < k; ++i) {
            d[v][i] = primal_vals[v * k + i];
        }
    }
    
    return {d, dual_obj};
}

std::tuple<std::vector<std::pair<int, int>>, double> calinescu_karloff_rabani_rounding(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& costs,
    const std::vector<int>& terminals,
    const std::vector<std::vector<double>>& d,
    int trials = 500
) {
    int k = terminals.size();
    double best_cost = std::numeric_limits<double>::infinity();
    std::vector<std::pair<int, int>> best_cut;
    
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, 0.5);
    
    for (int t = 0; t < trials; ++t) {
        std::vector<int> perm(k);
        for (int i = 0; i < k; ++i) perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), gen);
        
        double r = dist(gen);
        
        std::vector<int> assignment(n, -1);
        for (int v = 0; v < n; ++v) {
            bool assigned = false;
            for (int idx : perm) {
                if (d[v][idx] > r) {
                    assignment[v] = idx;
                    assigned = true;
                    break;
                }
            }
            if (!assigned) {
                assignment[v] = perm.back();
            }
        }
        
        std::vector<std::pair<int, int>> cut;
        double cost = 0.0;
        for (size_t e = 0; e < edges.size(); ++e) {
            int u = edges[e].first;
            int v = edges[e].second;
            if (assignment[u] != assignment[v]) {
                cut.push_back(edges[e]);
                cost += costs[e];
            }
        }
        
        if (cost < best_cost) {
            best_cost = cost;
            best_cut = cut;
        }
    }
    return {best_cut, best_cost};
}

} // namespace aal

using namespace aal;

void demo_multiway_cut_lp() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 19: Multiway Cut via LP Rounding\n");
    std::print("{:=^60}\n", "");
    
    int n = 6;
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
        {0, 5}, {2, 5}, {4, 5}
    };
    std::vector<double> costs = {2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0};
    std::vector<int> terminals = {0, 2, 4};
    
    std::print("\n1. Input Instance:\n");
    std::print("  Vertices: [0, 1, 2, 3, 4, 5]\n");
    std::print("  Terminals: [0, 2, 4]\n");
    std::print("  Edges & Costs:\n");
    for (size_t i = 0; i < edges.size(); ++i) {
        std::print("    ({}, {}) cost={:.1f}\n", edges[i].first, edges[i].second, costs[i]);
    }
    
    auto [d, lp_obj] = solve_multiway_cut_lp(n, edges, costs, terminals);
    std::print("\n  LP Relaxation Optimal Value: {:.4f}\n", lp_obj);
    std::print("  LP Vertex Embeddings on Simplex:\n");
    for (int v = 0; v < n; ++v) {
        std::print("    v_{}: [", v);
        for (size_t i = 0; i < d[v].size(); ++i) {
            if (i > 0) std::print(", ");
            std::print("{:.4f}", d[v][i]);
        }
        std::print("]\n");
    }
    
    auto [cut, cost] = calinescu_karloff_rabani_rounding(n, edges, costs, terminals, d, 500);
    std::print("\n2. CKR Rounding Results:\n");
    std::print("  Selected Cut Edges: [");
    for (size_t i = 0; i < cut.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("({}, {})", cut[i].first, cut[i].second);
    }
    std::print("]\n");
    std::print("  Cut Total Cost:     {:.2f}\n", cost);
}

#include "chapters.hpp"
#include "simplex.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <random>
#include <print>
#include <cmath>

namespace aal {

std::tuple<std::vector<double>, std::vector<double>, double> solve_max_sat_lp(
    int n_vars,
    const std::vector<std::pair<std::set<int>, std::set<int>>>& clauses,
    const std::vector<double>& weights
) {
    int n_clauses = clauses.size();
    int n_cols = n_vars + n_clauses;
    
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    
    // 1. z_j - sum_{i in C_j^+} y_i + sum_{i in C_j^-} y_i <= |C_j^-|
    for (int j = 0; j < n_clauses; ++j) {
        std::vector<double> row(n_cols, 0.0);
        row[n_vars + j] = 1.0;
        for (int i : clauses[j].first) {
            row[i] = -1.0;
        }
        for (int i : clauses[j].second) {
            row[i] = 1.0;
        }
        A.push_back(row);
        b.push_back(static_cast<double>(clauses[j].second.size()));
    }
    
    // 2. y_i <= 1
    for (int i = 0; i < n_vars; ++i) {
        std::vector<double> row(n_cols, 0.0);
        row[i] = 1.0;
        A.push_back(row);
        b.push_back(1.0);
    }
    
    // 3. z_j <= 1
    for (int j = 0; j < n_clauses; ++j) {
        std::vector<double> row(n_cols, 0.0);
        row[n_vars + j] = 1.0;
        A.push_back(row);
        b.push_back(1.0);
    }
    
    std::vector<double> c(n_cols, 0.0);
    for (int j = 0; j < n_clauses; ++j) {
        c[n_vars + j] = weights[j];
    }
    
    Simplex solver(A, b, c);
    auto [sol, obj] = solver.solve();
    
    if (sol.empty()) {
        return {std::vector<double>(n_vars, 0.0), std::vector<double>(n_clauses, 0.0), 0.0};
    }
    
    std::vector<double> y(sol.begin(), sol.begin() + n_vars);
    std::vector<double> z(sol.begin() + n_vars, sol.end());
    return {y, z, obj};
}

double evaluate_assignment(
    const std::vector<bool>& assignment,
    const std::vector<std::pair<std::set<int>, std::set<int>>>& clauses,
    const std::vector<double>& weights
) {
    double total_weight = 0.0;
    for (size_t j = 0; j < clauses.size(); ++j) {
        bool satisfied = false;
        for (int i : clauses[j].first) {
            if (assignment[i]) {
                satisfied = true;
                break;
            }
        }
        if (!satisfied) {
            for (int i : clauses[j].second) {
                if (!assignment[i]) {
                    satisfied = true;
                    break;
                }
            }
        }
        if (satisfied) {
            total_weight += weights[j];
        }
    }
    return total_weight;
}

std::tuple<std::vector<bool>, double, double> randomized_rounding_max_sat(
    int n_vars,
    const std::vector<std::pair<std::set<int>, std::set<int>>>& clauses,
    const std::vector<double>& weights,
    const std::vector<double>& y_lp,
    int trials = 200
) {
    double best_weight = -1.0;
    std::vector<bool> best_assignment;
    double total_weight = 0.0;
    
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int t = 0; t < trials; ++t) {
        std::vector<bool> assignment;
        for (int i = 0; i < n_vars; ++i) {
            assignment.push_back(dist(gen) < y_lp[i]);
        }
        
        double weight = evaluate_assignment(assignment, clauses, weights);
        total_weight += weight;
        if (weight > best_weight) {
            best_weight = weight;
            best_assignment = assignment;
        }
    }
    
    return {best_assignment, total_weight / trials, best_weight};
}

std::tuple<std::vector<bool>, double, double> coin_flip_max_sat(
    int n_vars,
    const std::vector<std::pair<std::set<int>, std::set<int>>>& clauses,
    const std::vector<double>& weights,
    int trials = 200
) {
    double best_weight = -1.0;
    std::vector<bool> best_assignment;
    double total_weight = 0.0;
    
    std::mt19937 gen(142);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int t = 0; t < trials; ++t) {
        std::vector<bool> assignment;
        for (int i = 0; i < n_vars; ++i) {
            assignment.push_back(dist(gen) < 0.5);
        }
        double weight = evaluate_assignment(assignment, clauses, weights);
        total_weight += weight;
        if (weight > best_weight) {
            best_weight = weight;
            best_assignment = assignment;
        }
    }
    
    return {best_assignment, total_weight / trials, best_weight};
}

} // namespace aal

using namespace aal;

void demo_randomized_rounding() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 16: Randomized Rounding for Max-SAT\n");
    std::print("{:=^60}\n", "");
    
    int n_vars = 3;
    std::vector<std::pair<std::set<int>, std::set<int>>> clauses = {
        {{0, 1}, {}},
        {{2}, {1}},
        {{}, {0, 2}},
        {{1}, {}}
    };
    std::vector<double> weights = {1.0, 2.0, 1.5, 3.0};
    
    std::print("\n1. Max-SAT CNF Instance:\n");
    std::print("  Variables: x_0, x_1, x_2\n");
    std::print("  Clauses & Weights:\n");
    for (size_t j = 0; j < clauses.size(); ++j) {
        std::string lits = "";
        bool first = true;
        for (int i : clauses[j].first) {
            if (!first) lits += " or ";
            lits += "x_" + std::to_string(i);
            first = false;
        }
        for (int i : clauses[j].second) {
            if (!first) lits += " or ";
            lits += "not x_" + std::to_string(i);
            first = false;
        }
        std::print("    C_{}: ({}), weight: {}\n", j, lits, weights[j]);
    }
    
    auto [y_lp, z_lp, lp_obj] = solve_max_sat_lp(n_vars, clauses, weights);
    std::print("\n  LP Optimal Objective Value: {:.4f}\n", lp_obj);
    std::print("  LP Variable Solution y*:    [");
    for (size_t i = 0; i < y_lp.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{:.4f}", y_lp[i]);
    }
    std::print("]\n");
    std::print("  LP Clause Solution z*:      [");
    for (size_t i = 0; i < z_lp.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{:.4f}", z_lp[i]);
    }
    std::print("]\n");
    
    auto [best_rr, avg_rr, max_rr] = randomized_rounding_max_sat(n_vars, clauses, weights, y_lp, 500);
    auto [best_cf, avg_cf, max_cf] = coin_flip_max_sat(n_vars, clauses, weights, 500);
    
    double opt_val = 0.0;
    std::vector<bool> best_all;
    for (int bitmask = 0; bitmask < 8; ++bitmask) {
        std::vector<bool> assignment;
        for (int i = 0; i < 3; ++i) {
            assignment.push_back((bitmask & (1 << i)) > 0);
        }
        double val = evaluate_assignment(assignment, clauses, weights);
        if (val > opt_val) {
            opt_val = val;
            best_all = assignment;
        }
    }
    
    std::print("\n  Exact Optimal Value:        {:.2f} (assignment: [", opt_val);
    for (size_t i = 0; i < best_all.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{}", best_all[i] ? "True" : "False");
    }
    std::print("])\n");
    std::print("  Randomized Rounding (avg):  {:.4f} (ratio to OPT: {:.4f}, bound: 1-1/e ~ 0.632)\n", avg_rr, avg_rr/opt_val);
    std::print("  Randomized Rounding (best): {:.4f}\n", max_rr);
    std::print("  Coin Flip Baseline (avg):   {:.4f} (ratio to OPT: {:.4f}, bound: 0.5)\n", avg_cf, avg_cf/opt_val);
    std::print("  Coin Flip Baseline (best):  {:.4f}\n", max_cf);
}

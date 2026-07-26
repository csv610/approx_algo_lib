#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <cmath>
#include <print>

namespace aal {

std::vector<bool> evaluate_sat_assignment(
    const std::vector<bool>& assignment,
    const std::vector<std::pair<int, int>>& clauses
) {
    std::vector<bool> satisfied;
    for (const auto& [l1, l2] : clauses) {
        bool sat1 = (l1 > 0) ? assignment[std::abs(l1)] : !assignment[std::abs(l1)];
        if (l2 == 0) {
            satisfied.push_back(sat1);
        } else {
            bool sat2 = (l2 > 0) ? assignment[std::abs(l2)] : !assignment[std::abs(l2)];
            satisfied.push_back(sat1 || sat2);
        }
    }
    return satisfied;
}

std::vector<std::vector<double>> optimize_max_2sat_vectors(
    int n_vars,
    const std::vector<std::pair<int, int>>& clauses,
    const std::vector<double>& weights,
    int dim = 8,
    double lr = 0.05,
    int epochs = 300
) {
    int n_vectors = n_vars + 1;
    
    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    
    std::vector<std::vector<double>> v(n_vectors, std::vector<double>(dim));
    for (int i = 0; i < n_vectors; ++i) {
        double mag2 = 0.0;
        for (int k = 0; k < dim; ++k) {
            v[i][k] = dist(gen);
            mag2 += v[i][k] * v[i][k];
        }
        double mag = std::sqrt(mag2);
        for (int k = 0; k < dim; ++k) v[i][k] /= mag;
    }
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        std::vector<std::vector<double>> new_v(n_vectors, std::vector<double>(dim));
        for (int i = 0; i < n_vectors; ++i) {
            std::vector<double> grad(dim, 0.0);
            
            for (size_t j = 0; j < clauses.size(); ++j) {
                int l1 = clauses[j].first;
                int l2 = clauses[j].second;
                double w = weights[j];
                
                int s1 = l1 > 0 ? 1 : -1;
                int v1_idx = std::abs(l1);
                
                if (l2 == 0) {
                    if (i == 0) {
                        for (int k = 0; k < dim; ++k) grad[k] += 0.5 * w * s1 * v[v1_idx][k];
                    } else if (i == v1_idx) {
                        for (int k = 0; k < dim; ++k) grad[k] += 0.5 * w * s1 * v[0][k];
                    }
                } else {
                    int s2 = l2 > 0 ? 1 : -1;
                    int v2_idx = std::abs(l2);
                    
                    if (i == 0) {
                        for (int k = 0; k < dim; ++k) grad[k] += 0.25 * w * (s1 * v[v1_idx][k] + s2 * v[v2_idx][k]);
                    } else if (i == v1_idx) {
                        for (int k = 0; k < dim; ++k) grad[k] += 0.25 * w * (s1 * v[0][k] - s1 * s2 * v[v2_idx][k]);
                    } else if (i == v2_idx) {
                        for (int k = 0; k < dim; ++k) grad[k] += 0.25 * w * (s2 * v[0][k] - s1 * s2 * v[v1_idx][k]);
                    }
                }
            }
            
            double mag2 = 0.0;
            std::vector<double> updated(dim);
            for (int k = 0; k < dim; ++k) {
                updated[k] = v[i][k] + lr * grad[k];
                mag2 += updated[k] * updated[k];
            }
            double mag = std::sqrt(mag2);
            if (mag < 1e-9) {
                new_v[i] = v[i];
            } else {
                for (int k = 0; k < dim; ++k) {
                    new_v[i][k] = updated[k] / mag;
                }
            }
        }
        v = new_v;
    }
    return v;
}

std::tuple<std::vector<bool>, double, double> goemans_williamson_max_2sat(
    int n_vars,
    const std::vector<std::pair<int, int>>& clauses,
    const std::vector<double>& weights,
    const std::vector<std::vector<double>>& vectors,
    int trials = 500
) {
    int dim = vectors[0].size();
    double best_weight = -1.0;
    std::vector<bool> best_assignment;
    double total_weight = 0.0;
    
    std::mt19937 gen(142);
    std::normal_distribution<double> dist(0.0, 1.0);
    
    for (int t = 0; t < trials; ++t) {
        std::vector<double> r(dim);
        double mag2 = 0.0;
        for (int k = 0; k < dim; ++k) {
            r[k] = dist(gen);
            mag2 += r[k] * r[k];
        }
        double mag = std::sqrt(mag2);
        for (int k = 0; k < dim; ++k) r[k] /= mag;
        
        double dot_v0 = 0.0;
        for (int k = 0; k < dim; ++k) dot_v0 += vectors[0][k] * r[k];
        bool sign_v0 = dot_v0 >= 0.0;
        
        std::vector<bool> assignment(n_vars + 1, false);
        for (int i = 1; i <= n_vars; ++i) {
            double dot_vi = 0.0;
            for (int k = 0; k < dim; ++k) dot_vi += vectors[i][k] * r[k];
            bool sign_vi = dot_vi >= 0.0;
            assignment[i] = (sign_vi == sign_v0);
        }
        
        auto sat_mask = evaluate_sat_assignment(assignment, clauses);
        double weight = 0.0;
        for (size_t j = 0; j < sat_mask.size(); ++j) {
            if (sat_mask[j]) weight += weights[j];
        }
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

void demo_max_2sat() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 26: SDP for Maximum 2-SAT\n");
    std::print("{:=^60}\n", "");
    
    int n_vars = 2;
    std::vector<std::pair<int, int>> clauses = {{1, 2}, {-1, 2}, {1, -2}, {-1, -2}};
    std::vector<double> weights = {1.0, 1.0, 1.0, 1.0};
    
    std::print("\n1. 2-Variable Max 2-SAT (All Combinations):\n");
    std::print("  Clauses & Weights:\n");
    for (size_t j = 0; j < clauses.size(); ++j) {
        int l1 = clauses[j].first;
        int l2 = clauses[j].second;
        std::string lit1 = l1 > 0 ? "x_" + std::to_string(l1) : "not x_" + std::to_string(-l1);
        std::string lit2 = l2 > 0 ? "x_" + std::to_string(l2) : "not x_" + std::to_string(-l2);
        std::print("    C_{}: ({} or {}), weight: {}\n", j, lit1, lit2, weights[j]);
    }
    
    auto vectors = optimize_max_2sat_vectors(n_vars, clauses, weights, 3, 0.1, 300);
    std::print("\n  Optimized Vector Embeddings:\n");
    std::print("    v_0 (True Anchor): [{:.4f}, {:.4f}, {:.4f}]\n", vectors[0][0], vectors[0][1], vectors[0][2]);
    for (int i = 1; i <= n_vars; ++i) {
        std::print("    v_{} (Var x_{}):   [{:.4f}, {:.4f}, {:.4f}]\n", i, i, vectors[i][0], vectors[i][1], vectors[i][2]);
    }
    
    auto [best_assign, avg_w, max_w] = goemans_williamson_max_2sat(n_vars, clauses, weights, vectors, 1000);
    std::print("\n  Exact Optimal Max 2-SAT Weight: 3.0\n");
    std::print("  GW Rounding Average Weight:    {:.4f} (ratio to OPT: {:.4f}, bound: 0.878)\n", avg_w, avg_w / 3.0);
    std::print("  GW Rounding Best Weight:       {:.4f} (assignment: [", max_w);
    for (size_t i = 1; i < best_assign.size(); ++i) {
        if (i > 1) std::print(", ");
        std::print("{}", best_assign[i] ? "True" : "False");
    }
    std::print("])\n");
    
    int n_vars2 = 5;
    std::vector<std::pair<int, int>> clauses2 = {
        {1, 2}, {-2, 3}, {-3, 4}, {4, 5}, {-5, -1},
        {1, 0}, {-2, 0}, {3, -5}
    };
    std::vector<double> weights2 = {1.5, 2.0, 1.0, 1.5, 2.0, 3.0, 1.0, 2.5};
    
    std::print("\n2. Larger Max 2-SAT Instance (5 variables, 8 clauses):\n");
    auto vectors2 = optimize_max_2sat_vectors(n_vars2, clauses2, weights2, 8, 0.05, 300);
    auto [best_assign2, avg_w2, max_w2] = goemans_williamson_max_2sat(n_vars2, clauses2, weights2, vectors2, 1000);
    
    double opt_w2 = 0.0;
    for (int bitmask = 0; bitmask < 32; ++bitmask) {
        std::vector<bool> assignment(6, false);
        for (int i = 0; i < 5; ++i) {
            assignment[i + 1] = (bitmask & (1 << i)) > 0;
        }
        auto sat_mask = evaluate_sat_assignment(assignment, clauses2);
        double w = 0.0;
        for (size_t j = 0; j < sat_mask.size(); ++j) {
            if (sat_mask[j]) w += weights2[j];
        }
        if (w > opt_w2) opt_w2 = w;
    }
    
    std::print("  Exact Optimal Max 2-SAT Weight: {:.2f}\n", opt_w2);
    std::print("  GW Rounding Average Weight:    {:.4f} (ratio to OPT: {:.4f}, bound: 0.878)\n", avg_w2, avg_w2 / opt_w2);
    std::print("  GW Rounding Best Weight:       {:.4f} (assignment: [", max_w2);
    for (size_t i = 1; i < best_assign2.size(); ++i) {
        if (i > 1) std::print(", ");
        std::print("{}", best_assign2[i] ? "True" : "False");
    }
    std::print("])\n");
}

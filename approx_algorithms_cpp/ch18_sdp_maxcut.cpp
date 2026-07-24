#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <cmath>
#include <print>

namespace aal {

double compute_cut_weight(
    const std::vector<int>& assignment,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& weights
) {
    double weight = 0.0;
    for (size_t i = 0; i < edges.size(); ++i) {
        int u = edges[i].first;
        int v = edges[i].second;
        if (assignment[u] != assignment[v]) {
            weight += weights[i];
        }
    }
    return weight;
}

std::vector<std::vector<double>> optimize_max_cut_vectors(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& weights,
    int dim = 8,
    double lr = 0.05,
    int epochs = 200
) {
    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    
    std::vector<std::vector<double>> v(n, std::vector<double>(dim));
    for (int i = 0; i < n; ++i) {
        double mag2 = 0.0;
        for (int k = 0; k < dim; ++k) {
            v[i][k] = dist(gen);
            mag2 += v[i][k] * v[i][k];
        }
        double mag = std::sqrt(mag2);
        for (int k = 0; k < dim; ++k) {
            v[i][k] /= mag;
        }
    }
    
    for (int epoch = 0; epoch < epochs; ++epoch) {
        std::vector<std::vector<double>> new_v(n, std::vector<double>(dim));
        for (int i = 0; i < n; ++i) {
            std::vector<double> grad(dim, 0.0);
            for (size_t e = 0; e < edges.size(); ++e) {
                int u = edges[e].first;
                int val_v = edges[e].second;
                double weight = weights[e];
                
                if (u == i) {
                    for (int k = 0; k < dim; ++k) grad[k] += weight * v[val_v][k];
                } else if (val_v == i) {
                    for (int k = 0; k < dim; ++k) grad[k] += weight * v[u][k];
                }
            }
            
            double mag2 = 0.0;
            std::vector<double> updated(dim);
            for (int k = 0; k < dim; ++k) {
                updated[k] = v[i][k] - lr * grad[k];
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

std::tuple<std::vector<int>, double, double> goemans_williamson_max_cut(
    int n,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& weights,
    const std::vector<std::vector<double>>& vectors,
    int trials = 500
) {
    int dim = vectors[0].size();
    double best_weight = -1.0;
    std::vector<int> best_assignment;
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
        
        std::vector<int> assignment(n);
        for (int i = 0; i < n; ++i) {
            double dot = 0.0;
            for (int k = 0; k < dim; ++k) {
                dot += vectors[i][k] * r[k];
            }
            assignment[i] = (dot >= 0.0) ? 1 : 0;
        }
        
        double weight = compute_cut_weight(assignment, edges, weights);
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

void demo_sdp_max_cut() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 18: Semidefinite Programming (Max-Cut)\n");
    std::print("{:=^60}\n", "");
    
    int n1 = 5;
    std::vector<std::pair<int, int>> edges1 = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}};
    std::vector<double> weights1 = {1.0, 1.0, 1.0, 1.0, 1.0};
    
    std::print("\n1. 5-Cycle Graph C_5 (unit weights):\n");
    std::print("  Vertices: [0, 1, 2, 3, 4]\n");
    std::print("  Edges:    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]\n");
    
    auto vectors1 = optimize_max_cut_vectors(n1, edges1, weights1, 2, 0.1, 300);
    
    std::print("\n  Optimized 2D Vector Embeddings on Circle:\n");
    for (int i = 0; i < n1; ++i) {
        std::print("    v_{}: [{:.4f}, {:.4f}]\n", i, vectors1[i][0], vectors1[i][1]);
    }
    
    auto [best_cut, avg_cut, max_cut] = goemans_williamson_max_cut(n1, edges1, weights1, vectors1, 1000);
    
    std::print("\n  Exact Optimal Max-Cut Value:  4.0\n");
    std::print("  GW Rounding Average Cut Value: {:.4f} (ratio to OPT: {:.4f}, bound: 0.878)\n", avg_cut, avg_cut/4.0);
    std::print("  GW Rounding Best Cut Value:    {:.4f} (assignment: [", max_cut);
    for (size_t i = 0; i < best_cut.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{}", best_cut[i]);
    }
    std::print("])\n");
    
    int n2 = 10;
    std::vector<std::pair<int, int>> edges2 = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
        {5, 7}, {7, 9}, {9, 6}, {6, 8}, {8, 5},
        {0, 5}, {1, 6}, {2, 7}, {3, 8}, {4, 9}
    };
    std::vector<double> weights2(edges2.size(), 1.0);
    
    std::print("\n2. Petersen Graph (10 vertices, 15 edges):\n");
    auto vectors2 = optimize_max_cut_vectors(n2, edges2, weights2, 8, 0.05, 300);
    auto [best_cut2, avg_cut2, max_cut2] = goemans_williamson_max_cut(n2, edges2, weights2, vectors2, 1000);
    
    std::print("  Exact Optimal Max-Cut Value:  12.0\n");
    std::print("  GW Rounding Average Cut Value: {:.4f} (ratio to OPT: {:.4f}, bound: 0.878)\n", avg_cut2, avg_cut2/12.0);
    std::print("  GW Rounding Best Cut Value:    {:.4f} (assignment: [", max_cut2);
    for (size_t i = 0; i < best_cut2.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{}", best_cut2[i]);
    }
    std::print("])\n");
}

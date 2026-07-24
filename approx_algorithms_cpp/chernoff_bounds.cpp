#include "chapters.hpp"
#include "simplex.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <random>
#include <cmath>
#include <algorithm>
#include <print>

namespace aal {

std::tuple<std::vector<double>, double> solve_set_cover_lp(
    const std::set<int>& universe,
    const std::map<int, std::set<int>>& sets,
    const std::map<int, double>& costs
) {
    std::vector<int> elements(universe.begin(), universe.end());
    std::sort(elements.begin(), elements.end());
    
    int n_sets = sets.size();
    int n_elements = elements.size();
    
    std::vector<std::vector<double>> A_dual;
    std::vector<double> b_dual;
    std::vector<double> c_dual(n_elements, 1.0);
    
    for (const auto& [s_idx, s_elements] : sets) {
        std::vector<double> row(n_elements, 0.0);
        for (int j = 0; j < n_elements; ++j) {
            if (s_elements.contains(elements[j])) {
                row[j] = 1.0;
            }
        }
        A_dual.push_back(row);
        b_dual.push_back(costs.at(s_idx));
    }
    
    Simplex solver(A_dual, b_dual, c_dual);
    auto [dual_sol, dual_obj] = solver.solve();
    
    if (dual_sol.empty()) {
        return {std::vector<double>(n_sets, 0.0), 0.0};
    }
    
    std::vector<double> primal_sol;
    for (int s_idx = 0; s_idx < n_sets; ++s_idx) {
        double val = solver.obj_row[n_elements + s_idx];
        primal_sol.push_back(std::max(0.0, val));
    }
    
    return {primal_sol, dual_obj};
}

std::tuple<std::set<int>, double, bool> set_cover_randomized_rounding(
    const std::set<int>& universe,
    const std::map<int, std::set<int>>& sets,
    const std::map<int, double>& costs,
    const std::vector<double>& x_lp,
    double c_factor = 1.5
) {
    int n_elements = universe.size();
    int n_sets = sets.size();
    
    int t = static_cast<int>(std::ceil(c_factor * std::log(std::max(2, n_elements))));
    
    std::set<int> chosen_sets;
    
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int i = 0; i < t; ++i) {
        for (int s_idx = 0; s_idx < n_sets; ++s_idx) {
            if (dist(gen) < x_lp[s_idx]) {
                chosen_sets.insert(s_idx);
            }
        }
    }
    
    std::set<int> covered;
    double cost = 0.0;
    for (int s_idx : chosen_sets) {
        for (int elem : sets.at(s_idx)) {
            covered.insert(elem);
        }
        cost += costs.at(s_idx);
    }
    
    bool is_valid = (covered == universe);
    return {chosen_sets, cost, is_valid};
}

} // namespace aal

using namespace aal;

void demo_chernoff_bounds() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 17: Set Cover via Randomized Rounding\n");
    std::print("{:=^60}\n", "");
    
    std::set<int> universe;
    for (int i = 0; i < 15; ++i) universe.insert(i);
    
    std::map<int, std::set<int>> sets = {
        {0, {0, 1, 2, 3, 4}},
        {1, {3, 4, 5, 6, 7}},
        {2, {6, 7, 8, 9, 10}},
        {3, {9, 10, 11, 12, 13}},
        {4, {12, 13, 14, 0, 1}},
        {5, {2, 5, 8, 11, 14}}
    };
    std::map<int, double> costs = {
        {0, 2.0}, {1, 2.0}, {2, 2.0}, {3, 2.0}, {4, 2.0}, {5, 1.5}
    };
    
    std::print("\n1. Set Cover Instance:\n");
    std::print("  Universe Size: {}\n", universe.size());
    std::print("  Available Sets & Costs:\n");
    for (const auto& [k, v] : sets) {
        std::print("    Set {}: [", k);
        bool first = true;
        for (int elem : v) {
            if (!first) std::print(", ");
            std::print("{}", elem);
            first = false;
        }
        std::print("] (cost: {:.1f})\n", costs.at(k));
    }
    
    auto [x_lp, lp_obj] = solve_set_cover_lp(universe, sets, costs);
    std::print("\n  LP Relaxation Optimal Value: {:.4f}\n", lp_obj);
    std::print("  LP Variable Solution x*:     [");
    for (size_t i = 0; i < x_lp.size(); ++i) {
        if (i > 0) std::print(", ");
        std::print("{:.4f}", x_lp[i]);
    }
    std::print("]\n");
    
    std::print("\n2. Randomized Rounding Simulation (varying scaling factor c):\n");
    std::vector<double> c_factors = {0.5, 1.0, 1.5, 2.0};
    for (double c : c_factors) {
        int t = static_cast<int>(std::ceil(c * std::log(universe.size())));
        int success_count = 0;
        double total_cost = 0.0;
        int n_trials = 200;
        
        for (int i = 0; i < n_trials; ++i) {
            auto [chosen, cost, is_valid] = set_cover_randomized_rounding(universe, sets, costs, x_lp, c);
            if (is_valid) {
                success_count++;
                total_cost += cost;
            }
        }
        
        double avg_cost = total_cost / std::max(1, success_count);
        double success_rate = static_cast<double>(success_count) / n_trials;
        std::print("  c={:.1f} (rounds={}): success rate={:.2f}%, avg cost of valid covers={:.2f}\n", c, t, success_rate * 100.0, avg_cost);
    }
}

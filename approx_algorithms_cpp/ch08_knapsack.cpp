#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iomanip>

namespace aal {

std::pair<std::vector<int>, int> knapsack_dp(const std::vector<int>& weights, const std::vector<int>& values, int capacity) {
    int n = weights.size();
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(capacity + 1, 0));
    
    for (int i = 1; i <= n; ++i) {
        for (int w = 0; w <= capacity; ++w) {
            if (weights[i - 1] <= w) {
                dp[i][w] = std::max(dp[i - 1][w], dp[i - 1][w - weights[i - 1]] + values[i - 1]);
            } else {
                dp[i][w] = dp[i - 1][w];
            }
        }
    }
    
    std::vector<int> selected;
    int w = capacity;
    for (int i = n; i > 0; --i) {
        if (dp[i][w] != dp[i - 1][w]) {
            selected.push_back(i - 1);
            w -= weights[i - 1];
        }
    }
    
    std::reverse(selected.begin(), selected.end());
    return {selected, dp[n][capacity]};
}

std::pair<std::vector<int>, int> knapsack_fptas(const std::vector<int>& weights, const std::vector<int>& values, int capacity, double epsilon) {
    int n = weights.size();
    if (n == 0) return {{}, 0};
    
    int v_max = *std::max_element(values.begin(), values.end());
    if (v_max == 0) return {{}, 0};
    
    double K = (epsilon * v_max) / n;
    if (K == 0) K = 1;
    
    std::vector<int> scaled_values(n);
    int max_scaled_val = 0;
    for (int i = 0; i < n; ++i) {
        scaled_values[i] = static_cast<int>(values[i] / K);
        max_scaled_val += scaled_values[i];
    }
    
    const int INF = 1e9;
    std::vector<int> dp(max_scaled_val + 1, INF);
    dp[0] = 0;
    
    std::vector<std::vector<bool>> choice(n + 1, std::vector<bool>(max_scaled_val + 1, false));
    
    for (int i = 1; i <= n; ++i) {
        int v = scaled_values[i - 1];
        int w = weights[i - 1];
        for (int val = max_scaled_val; val >= v; --val) {
            if (dp[val - v] != INF && dp[val - v] + w < dp[val]) {
                dp[val] = dp[val - v] + w;
                choice[i][val] = true;
            }
        }
    }
    
    int best_val = 0;
    for (int val = 0; val <= max_scaled_val; ++val) {
        if (dp[val] <= capacity) {
            best_val = val;
        }
    }
    
    std::vector<int> selected;
    int val = best_val;
    for (int i = n; i > 0; --i) {
        if (choice[i][val]) {
            selected.push_back(i - 1);
            val -= scaled_values[i - 1];
        }
    }
    
    int actual_value = 0;
    for (int idx : selected) {
        actual_value += values[idx];
    }
    
    std::reverse(selected.begin(), selected.end());
    return {selected, actual_value};
}

std::pair<std::vector<int>, int> knapsack_greedy_ratio(const std::vector<int>& weights, const std::vector<int>& values, int capacity) {
    int n = weights.size();
    std::vector<std::pair<int, double>> items(n);
    for (int i = 0; i < n; ++i) {
        double ratio = (weights[i] > 0) ? static_cast<double>(values[i]) / weights[i] : 1e9;
        items[i] = {i, ratio};
    }
    
    std::sort(items.begin(), items.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });
    
    std::vector<int> selected;
    int total_weight = 0;
    int total_value = 0;
    
    for (const auto& item : items) {
        int i = item.first;
        if (total_weight + weights[i] <= capacity) {
            selected.push_back(i);
            total_weight += weights[i];
            total_value += values[i];
        }
    }
    
    return {selected, total_value};
}

std::pair<std::vector<int>, int> knapsack_unbounded_dp(const std::vector<int>& weights, const std::vector<int>& values, int capacity) {
    int n = weights.size();
    std::vector<int> dp(capacity + 1, 0);
    std::vector<int> choice(capacity + 1, -1);
    
    for (int w = 1; w <= capacity; ++w) {
        for (int i = 0; i < n; ++i) {
            if (weights[i] <= w && dp[w - weights[i]] + values[i] > dp[w]) {
                dp[w] = dp[w - weights[i]] + values[i];
                choice[w] = i;
            }
        }
    }
    
    std::vector<int> selected;
    int w = capacity;
    while (w > 0 && choice[w] != -1) {
        selected.push_back(choice[w]);
        w -= weights[choice[w]];
    }
    
    return {selected, dp[capacity]};
}

} // namespace aal

using namespace aal;

void demo_knapsack() {
    std::cout << "============================================================\n";
    std::cout << "Chapter 8: Knapsack FPTAS\n";
    std::cout << "============================================================\n";
    
    std::vector<int> weights = {10, 20, 30, 40, 50};
    std::vector<int> values = {60, 100, 120, 200, 250};
    int capacity = 100;
    
    std::cout << "\nInstance: n=" << weights.size() << ", capacity=" << capacity << "\nWeights: [";
    for (size_t i = 0; i < weights.size(); ++i) std::cout << weights[i] << (i + 1 == weights.size() ? "" : ", ");
    std::cout << "]\nValues:  [";
    for (size_t i = 0; i < values.size(); ++i) std::cout << values[i] << (i + 1 == values.size() ? "" : ", ");
    std::cout << "]\n";
    
    auto [sel_exact, val_exact] = knapsack_dp(weights, values, capacity);
    std::cout << "\nExact DP: items=[";
    for (size_t i = 0; i < sel_exact.size(); ++i) std::cout << sel_exact[i] << (i + 1 == sel_exact.size() ? "" : ", ");
    std::cout << "], value=" << val_exact << "\n";
    
    std::vector<double> epsilons = {0.5, 0.25, 0.1, 0.05, 0.01};
    for (double eps : epsilons) {
        auto [sel, val] = knapsack_fptas(weights, values, capacity, eps);
        double ratio = val_exact > 0 ? static_cast<double>(val) / val_exact : 0.0;
        std::cout << "  FPTAS eps=" << std::fixed << std::setprecision(2) << eps 
                  << ": value=" << val << ", ratio=" << std::setprecision(4) << ratio 
                  << " (bound=" << (1.0 - eps) << ")\n";
    }
    
    std::cout << "\n--- Another Example ---\n";
    std::vector<int> weights2 = {2, 3, 4, 5};
    std::vector<int> values2 = {3, 4, 5, 6};
    int capacity2 = 8;
    
    auto [sel_exact2, val_exact2] = knapsack_dp(weights2, values2, capacity2);
    std::cout << "Exact: [";
    for (size_t i = 0; i < sel_exact2.size(); ++i) std::cout << sel_exact2[i] << (i + 1 == sel_exact2.size() ? "" : ", ");
    std::cout << "], value=" << val_exact2 << "\n";
    
    for (double eps : {0.2, 0.1}) {
        auto [sel, val] = knapsack_fptas(weights2, values2, capacity2, eps);
        std::cout << "  FPTAS eps=" << std::setprecision(1) << eps << ": [";
        for (size_t i = 0; i < sel.size(); ++i) std::cout << sel[i] << (i + 1 == sel.size() ? "" : ", ");
        std::cout << "], value=" << val << ", ratio=" << std::setprecision(4) << static_cast<double>(val)/val_exact2 << "\n";
    }
    
    std::cout << "\n--- Unbounded Knapsack ---\n";
    auto [sel_ub, val_ub] = knapsack_unbounded_dp(weights, values, capacity);
    std::cout << "Unbounded exact: items=[";
    for (size_t i = 0; i < sel_ub.size(); ++i) std::cout << sel_ub[i] << (i + 1 == sel_ub.size() ? "" : ", ");
    std::cout << "], value=" << val_ub << "\n";
}

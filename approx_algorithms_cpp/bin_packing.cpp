#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <map>

namespace aal {

std::vector<std::vector<double>> next_fit(const std::vector<double>& items, double capacity = 1.0) {
    std::vector<std::vector<double>> bins;
    std::vector<double> current_bin;
    double current_weight = 0.0;
    
    for (double item : items) {
        if (current_weight + item <= capacity) {
            current_bin.push_back(item);
            current_weight += item;
        } else {
            if (!current_bin.empty()) {
                bins.push_back(current_bin);
            }
            current_bin = {item};
            current_weight = item;
        }
    }
    
    if (!current_bin.empty()) {
        bins.push_back(current_bin);
    }
    return bins;
}

std::vector<std::vector<double>> first_fit(const std::vector<double>& items, double capacity = 1.0) {
    std::vector<std::vector<double>> bins;
    
    for (double item : items) {
        bool placed = false;
        for (auto& b : bins) {
            double sum = 0;
            for (double x : b) sum += x;
            if (sum + item <= capacity) {
                b.push_back(item);
                placed = true;
                break;
            }
        }
        if (!placed) {
            bins.push_back({item});
        }
    }
    
    return bins;
}

std::vector<std::vector<double>> first_fit_decreasing(const std::vector<double>& items, double capacity = 1.0) {
    std::vector<double> sorted_items = items;
    std::sort(sorted_items.begin(), sorted_items.end(), std::greater<double>());
    return first_fit(sorted_items, capacity);
}

std::vector<std::vector<int>> generate_configurations(const std::vector<double>& sizes, double cap = 1.0) {
    std::vector<std::vector<int>> configs;
    int n = sizes.size();
    
    auto backtrack = [&](auto& self, int idx, std::vector<int>& current_conf, double remaining_cap) -> void {
        if (idx == n) {
            int sum = 0;
            for (int x : current_conf) sum += x;
            if (sum > 0) {
                configs.push_back(current_conf);
            }
            return;
        }
        int max_count = static_cast<int>(remaining_cap / sizes[idx] + 1e-9);
        for (int count = 0; count <= max_count; ++count) {
            current_conf.push_back(count);
            self(self, idx + 1, current_conf, remaining_cap - count * sizes[idx]);
            current_conf.pop_back();
        }
    };
    
    std::vector<int> initial_conf;
    backtrack(backtrack, 0, initial_conf, cap);
    return configs;
}

std::vector<std::vector<double>> pack_large_dp(const std::vector<int>& counts, const std::vector<std::vector<int>>& configs, const std::vector<double>& sizes) {
    std::map<std::vector<int>, std::pair<int, std::vector<std::vector<double>>>> memo;
    
    auto solve = [&](auto& self, const std::vector<int>& state) -> std::pair<int, std::vector<std::vector<double>>> {
        int sum = 0;
        for (int x : state) sum += x;
        if (sum == 0) return {0, {}};
        
        if (memo.count(state)) return memo[state];
        
        int best_val = 1e9;
        std::vector<std::vector<double>> best_bins;
        
        for (const auto& conf : configs) {
            bool valid = true;
            for (size_t i = 0; i < state.size(); ++i) {
                if (state[i] < conf[i]) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                std::vector<int> next_state(state.size());
                for (size_t i = 0; i < state.size(); ++i) {
                    next_state[i] = state[i] - conf[i];
                }
                auto res = self(self, next_state);
                if (1 + res.first < best_val) {
                    best_val = 1 + res.first;
                    std::vector<double> new_bin;
                    for (size_t i = 0; i < conf.size(); ++i) {
                        for (int k = 0; k < conf[i]; ++k) {
                            new_bin.push_back(sizes[i]);
                        }
                    }
                    best_bins = {new_bin};
                    for (const auto& b : res.second) best_bins.push_back(b);
                }
            }
        }
        
        memo[state] = {best_val, best_bins};
        return memo[state];
    };
    
    return solve(solve, counts).second;
}

std::vector<std::vector<double>> bin_packing_aptas(std::vector<double> items, double eps = 0.3, double capacity = 1.0) {
    std::vector<double> large_items, small_items;
    for (double x : items) {
        if (x >= eps) large_items.push_back(x);
        else small_items.push_back(x);
    }
    
    if (large_items.empty()) {
        return first_fit(small_items, capacity);
    }
    
    std::sort(large_items.begin(), large_items.end());
    int n_large = large_items.size();
    int k = static_cast<int>(1.0 / (eps * eps));
    int q = n_large / k;
    
    std::vector<double> rounded_large;
    if (q == 0 || k == 0) {
        rounded_large = large_items;
    } else {
        std::vector<std::vector<double>> groups;
        for (int i = 0; i < k; ++i) {
            int start = i * q;
            int end = (i < k - 1) ? (i + 1) * q : n_large;
            std::vector<double> group(large_items.begin() + start, large_items.begin() + end);
            groups.push_back(group);
        }
        for (const auto& group : groups) {
            double max_size = *std::max_element(group.begin(), group.end());
            for (size_t i = 0; i < group.size(); ++i) {
                rounded_large.push_back(max_size);
            }
        }
    }
    
    std::vector<double> distinct_sizes;
    for (double val : rounded_large) {
        if (std::find(distinct_sizes.begin(), distinct_sizes.end(), val) == distinct_sizes.end()) {
            distinct_sizes.push_back(val);
        }
    }
    std::sort(distinct_sizes.begin(), distinct_sizes.end());
    
    std::vector<int> counts;
    for (double s : distinct_sizes) {
        counts.push_back(std::count(rounded_large.begin(), rounded_large.end(), s));
    }
    
    auto configs = generate_configurations(distinct_sizes, capacity);
    auto large_bins = pack_large_dp(counts, configs, distinct_sizes);
    
    std::vector<double> flattened_bins_items;
    for (const auto& b : large_bins) {
        for (double x : b) flattened_bins_items.push_back(x);
    }
    std::sort(flattened_bins_items.begin(), flattened_bins_items.end());
    
    std::map<double, std::vector<double>> item_map;
    for (size_t i = 0; i < flattened_bins_items.size(); ++i) {
        item_map[flattened_bins_items[i]].push_back(large_items[i]);
    }
    
    std::vector<std::vector<double>> final_large_bins;
    for (const auto& b : large_bins) {
        std::vector<double> new_bin;
        for (double item : b) {
            double orig_val = item_map[item].front();
            item_map[item].erase(item_map[item].begin());
            new_bin.push_back(orig_val);
        }
        final_large_bins.push_back(new_bin);
    }
    
    for (double item : small_items) {
        bool placed = false;
        for (auto& b : final_large_bins) {
            double sum = 0;
            for (double x : b) sum += x;
            if (sum + item <= capacity) {
                b.push_back(item);
                placed = true;
                break;
            }
        }
        if (!placed) {
            final_large_bins.push_back({item});
        }
    }
    
    return final_large_bins;
}

} // namespace aal

using namespace aal;

void print_bins(const std::vector<std::vector<double>>& bins) {
    std::cout << "[";
    for (size_t i = 0; i < bins.size(); ++i) {
        std::cout << "[";
        for (size_t j = 0; j < bins[i].size(); ++j) {
            std::cout << bins[i][j] << (j + 1 == bins[i].size() ? "" : ", ");
        }
        std::cout << "]" << (i + 1 == bins.size() ? "" : ", ");
    }
    std::cout << "]";
}

void demo_bin_packing() {
    std::cout << "============================================================\n";
    std::cout << "Chapter 9: Bin Packing Algorithms\n";
    std::cout << "============================================================\n";
    
    std::vector<double> items1 = {0.2, 0.5, 0.4, 0.7, 0.1, 0.3, 0.8};
    std::cout << "\n1. Input Items: [";
    for (size_t i = 0; i < items1.size(); ++i) std::cout << items1[i] << (i + 1 == items1.size() ? "" : ", ");
    std::cout << "]\n";
    
    auto nf1 = next_fit(items1);
    std::cout << "  Next-Fit (NF):           "; print_bins(nf1); std::cout << " (bins: " << nf1.size() << ")\n";
    
    auto ff1 = first_fit(items1);
    std::cout << "  First-Fit (FF):          "; print_bins(ff1); std::cout << " (bins: " << ff1.size() << ")\n";
    
    auto ffd1 = first_fit_decreasing(items1);
    std::cout << "  First-Fit Decreasing:    "; print_bins(ffd1); std::cout << " (bins: " << ffd1.size() << ")\n";
    
    auto aptas1 = bin_packing_aptas(items1, 0.4);
    std::cout << "  APTAS (eps=0.4):         "; print_bins(aptas1); std::cout << " (bins: " << aptas1.size() << ")\n";
    
    std::vector<double> items2 = {0.15, 0.33, 0.45, 0.12, 0.61, 0.38, 0.49, 0.52, 0.23, 0.29, 0.41, 0.19, 0.31, 0.27, 0.55, 0.48, 0.35, 0.11};
    std::cout << "\n2. Larger Instance (n=" << items2.size() << "):\n";
    std::cout << "  Next-Fit bins:           " << next_fit(items2).size() << "\n";
    std::cout << "  First-Fit bins:          " << first_fit(items2).size() << "\n";
    std::cout << "  First-Fit Decreasing:    " << first_fit_decreasing(items2).size() << "\n";
    std::cout << "  APTAS (eps=0.35) bins:   " << bin_packing_aptas(items2, 0.35).size() << "\n";
}

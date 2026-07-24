#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <map>

namespace aal {

std::vector<std::vector<double>> list_scheduling(const std::vector<double>& jobs, int m) {
    std::vector<std::vector<double>> schedule(m);
    std::vector<double> loads(m, 0.0);
    
    for (double job : jobs) {
        int min_idx = 0;
        for (int i = 1; i < m; ++i) {
            if (loads[i] < loads[min_idx]) {
                min_idx = i;
            }
        }
        schedule[min_idx].push_back(job);
        loads[min_idx] += job;
    }
    return schedule;
}

std::vector<std::vector<double>> lpt_scheduling(const std::vector<double>& jobs, int m) {
    std::vector<double> sorted_jobs = jobs;
    std::sort(sorted_jobs.begin(), sorted_jobs.end(), std::greater<double>());
    return list_scheduling(sorted_jobs, m);
}

// Ensure generate_configurations and pack_large_dp from ch09 are available.
// Redeclaring here inside aal to keep files independent (no external deps rule, stdlib only).
std::vector<std::vector<int>> generate_configurations_ms(const std::vector<double>& sizes, double cap) {
    std::vector<std::vector<int>> configs;
    int n = sizes.size();
    
    auto backtrack = [&](auto& self, int idx, std::vector<int>& current_conf, double remaining_cap) -> void {
        if (idx == n) {
            int sum = 0;
            for (int x : current_conf) sum += x;
            if (sum > 0) configs.push_back(current_conf);
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

std::vector<std::vector<double>> pack_large_dp_ms(const std::vector<int>& counts, const std::vector<std::vector<int>>& configs, const std::vector<double>& sizes) {
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
                if (state[i] < conf[i]) { valid = false; break; }
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
                        for (int k = 0; k < conf[i]; ++k) new_bin.push_back(sizes[i]);
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

std::pair<bool, std::vector<std::vector<double>>> check_schedule(const std::vector<double>& jobs, int m, double T, double eps) {
    std::vector<double> large_jobs, small_jobs;
    for (double p : jobs) {
        if (p > eps * T) large_jobs.push_back(p);
        else small_jobs.push_back(p);
    }
    
    if (large_jobs.empty()) {
        auto schedule = list_scheduling(small_jobs, m);
        double max_load = 0;
        for (const auto& sch : schedule) {
            double load = 0;
            for (double x : sch) load += x;
            max_load = std::max(max_load, load);
        }
        return {max_load <= (1 + eps) * T, schedule};
    }
    
    double delta = eps * eps * T;
    std::vector<double> rounded_large;
    std::vector<double> orig_large_sorted = large_jobs;
    std::sort(orig_large_sorted.begin(), orig_large_sorted.end());
    
    for (double p : orig_large_sorted) {
        double val = std::floor(p / delta) * delta;
        val = std::max(val, eps * T);
        rounded_large.push_back(val);
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
    
    auto configs = generate_configurations_ms(distinct_sizes, T);
    auto large_bins = pack_large_dp_ms(counts, configs, distinct_sizes);
    
    if (large_bins.size() > static_cast<size_t>(m)) {
        return {false, std::vector<std::vector<double>>(m)};
    }
    
    std::vector<double> flattened_bins_items;
    for (const auto& b : large_bins) {
        for (double x : b) flattened_bins_items.push_back(x);
    }
    std::sort(flattened_bins_items.begin(), flattened_bins_items.end());
    
    std::map<double, std::vector<double>> item_map;
    for (size_t i = 0; i < flattened_bins_items.size(); ++i) {
        item_map[flattened_bins_items[i]].push_back(orig_large_sorted[i]);
    }
    
    std::vector<std::vector<double>> schedule(m);
    for (size_t i = 0; i < large_bins.size(); ++i) {
        for (double item : large_bins[i]) {
            double orig_val = item_map[item].front();
            item_map[item].erase(item_map[item].begin());
            schedule[i].push_back(orig_val);
        }
    }
    
    std::vector<double> loads(m, 0.0);
    for (int i = 0; i < m; ++i) {
        for (double x : schedule[i]) loads[i] += x;
    }
    
    for (double job : small_jobs) {
        int min_idx = 0;
        for (int i = 1; i < m; ++i) {
            if (loads[i] < loads[min_idx]) min_idx = i;
        }
        schedule[min_idx].push_back(job);
        loads[min_idx] += job;
    }
    
    double max_load = 0;
    for (double l : loads) max_load = std::max(max_load, l);
    
    return {max_load <= (1 + eps) * T, schedule};
}

std::vector<std::vector<double>> makespan_ptas(const std::vector<double>& jobs, int m, double eps = 0.25) {
    double max_job = 0, sum_jobs = 0;
    for (double j : jobs) {
        max_job = std::max(max_job, j);
        sum_jobs += j;
    }
    double lb = std::max(max_job, sum_jobs / m);
    double ub = std::max(max_job, 2 * sum_jobs / m);
    
    std::vector<std::vector<double>> best_schedule;
    for (int step = 0; step < 15; ++step) {
        double mid = (lb + ub) / 2;
        auto res = check_schedule(jobs, m, mid, eps);
        if (res.first) {
            best_schedule = res.second;
            ub = mid;
        } else {
            lb = mid;
        }
    }
    
    if (best_schedule.empty()) {
        best_schedule = check_schedule(jobs, m, ub, eps).second;
    }
    return best_schedule;
}

} // namespace aal

using namespace aal;

double calc_makespan(const std::vector<std::vector<double>>& schedule) {
    double ms = 0;
    for (const auto& sch : schedule) {
        double sum = 0;
        for (double x : sch) sum += x;
        ms = std::max(ms, sum);
    }
    return ms;
}

void print_schedule(const std::vector<std::vector<double>>& sch) {
    std::cout << "[";
    for (size_t i = 0; i < sch.size(); ++i) {
        std::cout << "[";
        for (size_t j = 0; j < sch[i].size(); ++j) {
            std::cout << sch[i][j] << (j + 1 == sch[i].size() ? "" : ", ");
        }
        std::cout << "]" << (i + 1 == sch.size() ? "" : ", ");
    }
    std::cout << "]";
}

void demo_makespan() {
    std::cout << "============================================================\n";
    std::cout << "Chapter 10: Minimum Makespan Scheduling\n";
    std::cout << "============================================================\n";
    
    std::vector<double> jobs1 = {2.0, 3.0, 4.0, 6.0, 2.0};
    int m1 = 2;
    std::cout << "\n1. Jobs: [";
    for (size_t i = 0; i < jobs1.size(); ++i) std::cout << jobs1[i] << (i + 1 == jobs1.size() ? "" : ", ");
    std::cout << "] on " << m1 << " machines\n";
    
    auto sched_list = list_scheduling(jobs1, m1);
    std::cout << "  List Scheduling:      "; print_schedule(sched_list);
    std::cout << " (makespan: " << calc_makespan(sched_list) << ")\n";
    
    auto sched_lpt = lpt_scheduling(jobs1, m1);
    std::cout << "  LPT Heuristic:        "; print_schedule(sched_lpt);
    std::cout << " (makespan: " << calc_makespan(sched_lpt) << ")\n";
    
    auto sched_ptas = makespan_ptas(jobs1, m1, 0.25);
    std::cout << "  PTAS (eps=0.25):      "; print_schedule(sched_ptas);
    std::cout << " (makespan: " << calc_makespan(sched_ptas) << ")\n";
    
    std::vector<double> jobs2 = {1.2, 2.5, 3.1, 4.0, 1.8, 2.2, 5.0, 3.5, 0.9, 1.6};
    int m2 = 3;
    std::cout << "\n2. Larger Instance (n=" << jobs2.size() << ") on " << m2 << " machines:\n";
    
    std::cout << "  List Scheduling makespan:   " << calc_makespan(list_scheduling(jobs2, m2)) << "\n";
    std::cout << "  LPT Heuristic makespan:     " << calc_makespan(lpt_scheduling(jobs2, m2)) << "\n";
    std::cout << "  PTAS (eps=0.2) makespan:    " << calc_makespan(makespan_ptas(jobs2, m2, 0.2)) << "\n";
}

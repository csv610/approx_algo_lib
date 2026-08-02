/**
 * @file unrelated_scheduling.cpp
 * @brief Scheduling on Unrelated Parallel Machines (Vazirani, Ch. 17)
 *
 * Algorithms:
 * - Greedy: assign each job to its fastest machine
 * - Local search: swap-based improvement
 * - LP relaxation approximation (greedy-based)
 */

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <numeric>
#include "json.hpp"

using json = nlohmann::json;
using namespace std;

namespace aal {

using Graph = map<int, map<int, double>>;

struct SchedulingResult {
    double makespan;
    vector<int> assignment;
};

SchedulingResult greedy_schedule(const vector<vector<double>>& processing_times) {
    int m = processing_times.size();
    int n = processing_times[0].size();

    vector<int> assignment(n, 0);
    vector<double> machine_load(m, 0.0);

    for (int j = 0; j < n; ++j) {
        int best_machine = 0;
        double best_time = processing_times[0][j];
        for (int i = 1; i < m; ++i) {
            if (processing_times[i][j] < best_time) {
                best_time = processing_times[i][j];
                best_machine = i;
            }
        }
        assignment[j] = best_machine;
        machine_load[best_machine] += best_time;
    }

    double makespan = 0.0;
    for (double load : machine_load) {
        makespan = max(makespan, load);
    }

    return {makespan, assignment};
}

double compute_makespan(const vector<vector<double>>& processing_times,
                        const vector<int>& assignment) {
    int m = processing_times.size();
    vector<double> machine_load(m, 0.0);

    for (size_t j = 0; j < assignment.size(); ++j) {
        machine_load[assignment[j]] += processing_times[assignment[j]][j];
    }

    double makespan = 0.0;
    for (double load : machine_load) {
        makespan = max(makespan, load);
    }
    return makespan;
}

SchedulingResult local_search_schedule(const vector<vector<double>>& processing_times,
                                       const vector<int>& initial_assignment,
                                       int max_iterations = 1000) {
    int m = processing_times.size();
    int n = processing_times[0].size();

    vector<int> assignment = initial_assignment;
    vector<double> machine_load(m, 0.0);

    for (int j = 0; j < n; ++j) {
        machine_load[assignment[j]] += processing_times[assignment[j]][j];
    }

    bool improved = true;
    int iterations = 0;

    while (improved && iterations < max_iterations) {
        improved = false;
        ++iterations;

        for (int j1 = 0; j1 < n; ++j1) {
            for (int j2 = j1 + 1; j2 < n; ++j2) {
                int m1 = assignment[j1];
                int m2 = assignment[j2];

                if (m1 == m2) continue;

                double current_makespan = *max_element(machine_load.begin(), machine_load.end());

                machine_load[m1] -= processing_times[m1][j1];
                machine_load[m2] -= processing_times[m2][j2];

                machine_load[m1] += processing_times[m1][j2];
                machine_load[m2] += processing_times[m2][j1];

                double new_makespan = *max_element(machine_load.begin(), machine_load.end());

                if (new_makespan < current_makespan) {
                    swap(assignment[j1], assignment[j2]);
                    improved = true;
                } else {
                    machine_load[m1] -= processing_times[m1][j2];
                    machine_load[m2] -= processing_times[m2][j1];

                    machine_load[m1] += processing_times[m1][j1];
                    machine_load[m2] += processing_times[m2][j2];
                }
            }
        }
    }

    double makespan = *max_element(machine_load.begin(), machine_load.end());
    return {makespan, assignment};
}

SchedulingResult lp_relaxation_schedule(const vector<vector<double>>& processing_times) {
    int m = processing_times.size();
    int n = processing_times[0].size();

    vector<int> assignment(n, 0);
    vector<double> machine_load(m, 0.0);

    vector<double> lower_bound_candidates;

    for (int j = 0; j < n; ++j) {
        double min_time = processing_times[0][j];
        for (int i = 1; i < m; ++i) {
            min_time = min(min_time, processing_times[i][j]);
        }
        lower_bound_candidates.push_back(min_time);
    }

    double total_min = 0.0;
    for (double t : lower_bound_candidates) total_min += t;
    double lb = total_min / m;

    double max_single = 0.0;
    for (int j = 0; j < n; ++j) {
        double min_time = *min_element(processing_times.begin(), processing_times.end(),
                                       [j](const vector<double>& a, const vector<double>& b) {
                                           return a[j] < b[j];
                                       })[j];
        max_single = max(max_single, min_time);
    }
    lb = max(lb, max_single);

    auto greedy = greedy_schedule(processing_times);

    vector<int> sorted_jobs(n);
    iota(sorted_jobs.begin(), sorted_jobs.end(), 0);
    sort(sorted_jobs.begin(), sorted_jobs.end(), [&](int a, int b) {
        double min_a = processing_times[0][a];
        double min_b = processing_times[0][b];
        for (int i = 1; i < m; ++i) {
            min_a = min(min_a, processing_times[i][a]);
            min_b = min(min_b, processing_times[i][b]);
        }
        return min_a > min_b;
    });

    fill(machine_load.begin(), machine_load.end(), 0.0);
    for (int j : sorted_jobs) {
        int best_machine = 0;
        double best_load_time = machine_load[0] + processing_times[0][j];
        for (int i = 1; i < m; ++i) {
            double candidate = machine_load[i] + processing_times[i][j];
            if (candidate < best_load_time) {
                best_load_time = candidate;
                best_machine = i;
            }
        }
        assignment[j] = best_machine;
        machine_load[best_machine] += processing_times[best_machine][j];
    }

    double makespan = *max_element(machine_load.begin(), machine_load.end());

    if (makespan > greedy.makespan) {
        return greedy;
    }

    return {makespan, assignment};
}

void solve(istream& in, ostream& out) {
    json input;
    in >> input;

    vector<vector<double>> processing_times = input["processing_times"].get<vector<vector<double>>>();

    auto result_greedy = greedy_schedule(processing_times);
    auto result_ls = local_search_schedule(processing_times, result_greedy.assignment);
    auto result_lp = lp_relaxation_schedule(processing_times);

    json output;
    output["algorithm"] = "unrelated_scheduling";
    output["makespan_greedy"] = result_greedy.makespan;
    output["makespan_local_search"] = result_ls.makespan;
    output["makespan_lp_approx"] = result_lp.makespan;
    output["assignment_greedy"] = result_greedy.assignment;
    output["assignment_local_search"] = result_ls.assignment;
    output["assignment_lp_approx"] = result_lp.assignment;

    out << output.dump(2) << endl;
}

} // namespace aal

int main() {
    aal::solve(cin, cout);
    return 0;
}

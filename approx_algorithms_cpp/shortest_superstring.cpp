#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>

namespace aal {

int compute_overlap(const std::string& s1, const std::string& s2) {
    int max_len = std::min(s1.length(), s2.length());
    for (int l = max_len; l > 0; --l) {
        if (s1.compare(s1.length() - l, l, s2, 0, l) == 0) {
            return l;
        }
    }
    return 0;
}

std::vector<std::string> preprocess_substrings(const std::vector<std::string>& strings) {
    std::vector<std::string> sorted_strings = strings;
    std::sort(sorted_strings.begin(), sorted_strings.end(), [](const std::string& a, const std::string& b) {
        return a.length() > b.length();
    });
    std::vector<std::string> filtered;
    for (const auto& s : sorted_strings) {
        bool is_sub = false;
        for (const auto& other : filtered) {
            if (other.find(s) != std::string::npos) {
                is_sub = true;
                break;
            }
        }
        if (!is_sub) filtered.push_back(s);
    }
    return filtered;
}

std::pair<std::vector<int>, double> find_minimum_cycle_cover(const std::vector<std::vector<double>>& cost_matrix) {
    int n = cost_matrix.size();
    double best_cost = std::numeric_limits<double>::infinity();
    std::vector<int> best_perm;
    
    std::vector<int> perm;
    std::vector<bool> visited(n, false);
    
    auto backtrack = [&](auto& self, int curr, double cost) -> void {
        if (cost >= best_cost) return;
        if (curr == n) {
            best_cost = cost;
            best_perm = perm;
            return;
        }
        for (int next_val = 0; next_val < n; ++next_val) {
            if (!visited[next_val]) {
                visited[next_val] = true;
                perm.push_back(next_val);
                self(self, curr + 1, cost + cost_matrix[curr][next_val]);
                perm.pop_back();
                visited[next_val] = false;
            }
        }
    };
    
    backtrack(backtrack, 0, 0.0);
    return {best_perm, best_cost};
}

std::vector<std::vector<int>> extract_cycles(const std::vector<int>& perm) {
    int n = perm.size();
    std::vector<bool> visited(n, false);
    std::vector<std::vector<int>> cycles;
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            std::vector<int> cycle;
            int curr = i;
            while (!visited[curr]) {
                visited[curr] = true;
                cycle.push_back(curr);
                curr = perm[curr];
            }
            cycles.push_back(cycle);
        }
    }
    return cycles;
}

std::string greedy_superstring(const std::vector<std::string>& strings) {
    std::vector<std::string> T = strings;
    while (T.size() > 1) {
        int max_ov = -1;
        int best_i = -1, best_j = -1;
        for (size_t i = 0; i < T.size(); ++i) {
            for (size_t j = 0; j < T.size(); ++j) {
                if (i != j) {
                    int ov = compute_overlap(T[i], T[j]);
                    if (ov > max_ov) {
                        max_ov = ov;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }
        
        std::string s_i = T[best_i];
        std::string s_j = T[best_j];
        std::string merged = s_i + s_j.substr(max_ov);
        
        if (best_i < best_j) {
            T.erase(T.begin() + best_j);
            T.erase(T.begin() + best_i);
        } else {
            T.erase(T.begin() + best_i);
            T.erase(T.begin() + best_j);
        }
        T.push_back(merged);
    }
    return T.empty() ? "" : T[0];
}

std::vector<std::string> cycle_cover_to_strings(const std::vector<std::vector<int>>& cycles, const std::vector<std::string>& strings) {
    std::vector<std::string> cycle_strings;
    for (const auto& cycle : cycles) {
        if (cycle.size() == 1) {
            cycle_strings.push_back(strings[cycle[0]]);
            continue;
        }
        std::string merged = strings[cycle[0]];
        for (size_t idx = 0; idx < cycle.size() - 1; ++idx) {
            int u = cycle[idx];
            int v = cycle[idx + 1];
            int ov = compute_overlap(strings[u], strings[v]);
            merged += strings[v].substr(ov);
        }
        int last = cycle.back();
        int first = cycle.front();
        int ov = compute_overlap(strings[last], strings[first]);
        cycle_strings.push_back(merged);
    }
    return cycle_strings;
}

std::string shortest_superstring_4approx(const std::vector<std::string>& strings) {
    auto cleaned = preprocess_substrings(strings);
    if (cleaned.empty()) return "";
    if (cleaned.size() == 1) return cleaned[0];
    
    int n = cleaned.size();
    std::vector<std::vector<double>> cost_matrix(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cost_matrix[i][j] = cleaned[j].length() - compute_overlap(cleaned[i], cleaned[j]);
        }
    }
    
    auto [perm, _] = find_minimum_cycle_cover(cost_matrix);
    auto cycles = extract_cycles(perm);
    auto cycle_strings = cycle_cover_to_strings(cycles, cleaned);
    
    std::string result = "";
    for (const auto& s : cycle_strings) result += s;
    return result;
}

std::string shortest_superstring_3approx(const std::vector<std::string>& strings) {
    auto cleaned = preprocess_substrings(strings);
    if (cleaned.empty()) return "";
    if (cleaned.size() == 1) return cleaned[0];
    
    int n = cleaned.size();
    std::vector<std::vector<double>> cost_matrix(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cost_matrix[i][j] = cleaned[j].length() - compute_overlap(cleaned[i], cleaned[j]);
        }
    }
    
    auto [perm, _] = find_minimum_cycle_cover(cost_matrix);
    auto cycles = extract_cycles(perm);
    auto cycle_strings = cycle_cover_to_strings(cycles, cleaned);
    
    return greedy_superstring(cycle_strings);
}

} // namespace aal

void demo_shortest_superstring() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 7: Shortest Common Superstring\n";
    std::cout << "============================================================\n";
    
    std::vector<std::string> S1 = {"abc", "bcd", "cde", "def"};
    std::cout << "\n1. Input Strings: {abc, bcd, cde, def}\n";
    std::cout << "  Greedy Superstring:  " << greedy_superstring(S1) << "\n";
    std::cout << "  4-Approx Superstring: " << shortest_superstring_4approx(S1) << "\n";
    std::cout << "  3-Approx Superstring: " << shortest_superstring_3approx(S1) << "\n";
    
    std::vector<std::string> S2 = {"CATG", "ATGT", "TGTA", "GTAC", "TACA"};
    std::cout << "\n2. DNA Fragments: {CATG, ATGT, TGTA, GTAC, TACA}\n";
    std::cout << "  Greedy Superstring:  " << greedy_superstring(S2) << "\n";
    std::cout << "  4-Approx Superstring: " << shortest_superstring_4approx(S2) << "\n";
    std::cout << "  3-Approx Superstring: " << shortest_superstring_3approx(S2) << "\n";
}

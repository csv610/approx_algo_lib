#include "chapters.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <print>
#include <limits>
#include <random>

namespace aal {

struct Facility {
    double cost;
    std::map<int, double> clients;
};

std::tuple<std::set<int>, std::map<int, int>, double> facility_location_greedy(
    const std::map<int, Facility>& facilities,
    const std::vector<int>& clients
) {
    std::set<int> open_facilities;
    std::map<int, int> assignments;
    std::set<int> unassigned(clients.begin(), clients.end());
    double total_cost = 0.0;
    
    while (!unassigned.empty()) {
        int best_facility = -1;
        std::set<int> best_clients;
        double best_ratio = std::numeric_limits<double>::infinity();
        
        for (const auto& [i, fac] : facilities) {
            if (open_facilities.contains(i)) continue;
            
            std::vector<std::pair<int, double>> new_clients;
            for (int j : unassigned) {
                if (fac.clients.contains(j)) {
                    new_clients.push_back({j, fac.clients.at(j)});
                }
            }
            
            if (new_clients.empty()) continue;
            
            double conn_cost = 0.0;
            for (auto& p : new_clients) conn_cost += p.second;
            
            double total = fac.cost + conn_cost;
            double ratio = total / new_clients.size();
            
            if (ratio < best_ratio) {
                best_ratio = ratio;
                best_facility = i;
                best_clients.clear();
                for (auto& p : new_clients) best_clients.insert(p.first);
            }
        }
        
        if (best_facility == -1) {
            for (int j : unassigned) {
                double best_dist = std::numeric_limits<double>::infinity();
                int best_f = -1;
                for (int i : open_facilities) {
                    if (facilities.at(i).clients.contains(j)) {
                        double d = facilities.at(i).clients.at(j);
                        if (d < best_dist) {
                            best_dist = d;
                            best_f = i;
                        }
                    }
                }
                if (best_f != -1) {
                    assignments[j] = best_f;
                    total_cost += best_dist;
                }
            }
            break;
        }
        
        open_facilities.insert(best_facility);
        total_cost += facilities.at(best_facility).cost;
        
        for (int j : best_clients) {
            if (unassigned.contains(j)) {
                assignments[j] = best_facility;
                total_cost += facilities.at(best_facility).clients.at(j);
                unassigned.erase(j);
            }
        }
    }
    
    return {open_facilities, assignments, total_cost};
}

std::tuple<std::set<int>, std::map<int, int>, double> k_median_lp_rounding(
    const std::map<int, Facility>& facilities,
    const std::vector<int>& clients,
    int k
) {
    std::set<int> open_facilities;
    int count = 0;
    for (const auto& [i, fac] : facilities) {
        if (count < k) {
            open_facilities.insert(i);
            count++;
        }
    }
    
    auto compute_cost = [&](const std::set<int>& open_set) -> std::pair<std::map<int, int>, double> {
        std::map<int, int> assign;
        double total = 0.0;
        for (int j : clients) {
            double best_dist = std::numeric_limits<double>::infinity();
            int best_f = -1;
            for (int i : open_set) {
                double d = facilities.at(i).clients.contains(j) ? facilities.at(i).clients.at(j) : std::numeric_limits<double>::infinity();
                if (d < best_dist) {
                    best_dist = d;
                    best_f = i;
                }
            }
            if (best_f != -1) {
                assign[j] = best_f;
                total += best_dist;
            }
        }
        return {assign, total};
    };
    
    auto [assignments, cost] = compute_cost(open_facilities);
    
    bool improved = true;
    while (improved) {
        improved = false;
        std::set<int> closed;
        for (const auto& [i, fac] : facilities) {
            if (!open_facilities.contains(i)) {
                closed.insert(i);
            }
        }
        
        for (int o : open_facilities) {
            for (int c : closed) {
                std::set<int> new_open = open_facilities;
                new_open.erase(o);
                new_open.insert(c);
                auto [new_assign, new_cost] = compute_cost(new_open);
                if (new_cost < cost) {
                    open_facilities = new_open;
                    assignments = new_assign;
                    cost = new_cost;
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
    }
    
    return {open_facilities, assignments, cost};
}

std::tuple<std::set<int>, std::map<int, int>, double> facility_location_primal_dual(
    const std::map<int, Facility>& facilities,
    const std::vector<int>& clients
) {
    std::map<int, double> alpha;
    for (int j : clients) alpha[j] = 0.0;
    
    std::set<int> open_facilities;
    std::map<int, int> assignments;
    std::set<int> connected;
    
    std::set<int> active_facilities;
    for (const auto& [i, fac] : facilities) {
        active_facilities.insert(i);
    }
    
    while (connected.size() < clients.size()) {
        double min_delta = std::numeric_limits<double>::infinity();
        
        for (int i : active_facilities) {
            double sum_alpha = 0.0;
            for (int j : clients) {
                if (!connected.contains(j) && facilities.at(i).clients.contains(j)) {
                    sum_alpha += alpha[j];
                }
            }
            double slack = facilities.at(i).cost - sum_alpha;
            if (slack < min_delta) {
                min_delta = slack;
            }
        }
        
        for (int j : clients) {
            if (connected.contains(j)) continue;
            for (int i : active_facilities) {
                if (facilities.at(i).clients.contains(j)) {
                    double slack = facilities.at(i).clients.at(j) - alpha[j];
                    if (slack < min_delta) {
                        min_delta = slack;
                    }
                }
            }
        }
        
        if (min_delta == std::numeric_limits<double>::infinity() || min_delta <= 0) {
            break;
        }
        
        for (int j : clients) {
            if (!connected.contains(j)) {
                alpha[j] += min_delta;
            }
        }
        
        std::vector<int> act_fac(active_facilities.begin(), active_facilities.end());
        for (int i : act_fac) {
            double sum_alpha = 0.0;
            for (int j : clients) {
                if (!connected.contains(j) && facilities.at(i).clients.contains(j)) {
                    sum_alpha += alpha[j];
                }
            }
            if (sum_alpha >= facilities.at(i).cost - 1e-9) {
                open_facilities.insert(i);
                active_facilities.erase(i);
                for (int j : clients) {
                    if (!connected.contains(j) && facilities.at(i).clients.contains(j)) {
                        assignments[j] = i;
                        connected.insert(j);
                    }
                }
            }
        }
    }
    
    for (int j : clients) {
        if (!assignments.contains(j)) {
            double best_dist = std::numeric_limits<double>::infinity();
            int best_f = -1;
            for (int i : open_facilities) {
                if (facilities.at(i).clients.contains(j)) {
                    double d = facilities.at(i).clients.at(j);
                    if (d < best_dist) {
                        best_dist = d;
                        best_f = i;
                    }
                }
            }
            if (best_f != -1) {
                assignments[j] = best_f;
            }
        }
    }
    
    double total_cost = 0.0;
    for (int i : open_facilities) {
        total_cost += facilities.at(i).cost;
    }
    for (const auto& [j, i] : assignments) {
        total_cost += facilities.at(i).clients.at(j);
    }
    
    return {open_facilities, assignments, total_cost};
}

} // namespace aal

using namespace aal;

void demo_facility_location() {
    std::print("{:=^60}\n", "");
    std::print("Chapter 24: Facility Location\n");
    std::print("{:=^60}\n", "");
    
    std::map<int, Facility> facilities = {
        {0, {10.0, {{0, 2.0}, {1, 5.0}, {2, 3.0}, {3, 8.0}, {4, 7.0}}}},
        {1, {8.0,  {{0, 6.0}, {1, 2.0}, {2, 4.0}, {3, 3.0}, {4, 5.0}}}},
        {2, {12.0, {{0, 4.0}, {1, 6.0}, {2, 2.0}, {3, 5.0}, {4, 3.0}}}}
    };
    std::vector<int> clients = {0, 1, 2, 3, 4};
    
    std::print("\n1. Greedy Facility Location\n");
    auto [open_f, assign, cost] = facility_location_greedy(facilities, clients);
    std::print("  Open facilities: {{");
    bool first = true;
    for (int i : open_f) { if (!first) std::print(", "); std::print("{}", i); first = false; }
    std::print("}}\n  Assignments: {{");
    first = true;
    for (const auto& [k, v] : assign) { if (!first) std::print(", "); std::print("{}: {}", k, v); first = false; }
    std::print("}}\n  Total cost: {}\n", cost);
    
    std::print("\n2. Primal-Dual 3-approx\n");
    auto [open_f2, assign2, cost2] = facility_location_primal_dual(facilities, clients);
    std::print("  Open facilities: {{");
    first = true;
    for (int i : open_f2) { if (!first) std::print(", "); std::print("{}", i); first = false; }
    std::print("}}\n  Assignments: {{");
    first = true;
    for (const auto& [k, v] : assign2) { if (!first) std::print(", "); std::print("{}: {}", k, v); first = false; }
    std::print("}}\n  Total cost: {}\n", cost2);
    
    std::print("\n3. k-Median (k=2) Local Search\n");
    auto [open_f3, assign3, cost3] = k_median_lp_rounding(facilities, clients, 2);
    std::print("  Open facilities: {{");
    first = true;
    for (int i : open_f3) { if (!first) std::print(", "); std::print("{}", i); first = false; }
    std::print("}}\n  Assignments: {{");
    first = true;
    for (const auto& [k, v] : assign3) { if (!first) std::print(", "); std::print("{}: {}", k, v); first = false; }
    std::print("}}\n  Total cost: {}\n", cost3);
    
    std::print("\n--- Larger Example ---\n");
    std::mt19937 gen(42);
    std::uniform_int_distribution<> cost_dist(5, 20);
    std::uniform_int_distribution<> cli_dist(1, 10);
    
    int n_fac = 10;
    int n_cli = 20;
    std::map<int, Facility> facilities_large;
    for (int i = 0; i < n_fac; ++i) {
        Facility f;
        f.cost = cost_dist(gen);
        for (int j = 0; j < n_cli; ++j) {
            f.clients[j] = cli_dist(gen);
        }
        facilities_large[i] = f;
    }
    std::vector<int> clients_large;
    for (int j = 0; j < n_cli; ++j) clients_large.push_back(j);
    
    auto [open_f4, assign4, cost4] = facility_location_greedy(facilities_large, clients_large);
    std::print("  n_fac={}, n_cli={}\n", n_fac, n_cli);
    std::print("  Open: {} facilities, Cost: {:.1f}\n", open_f4.size(), cost4);
}

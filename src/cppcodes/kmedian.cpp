// =============================================================================
// k-Median Problem (Vazirani Chapter 25)
// Algorithms: Greedy, 1-Swap Local Search, Simple Randomized
// =============================================================================
// Given a set of clients and potential facility locations with distances,
// open exactly k facilities to minimize the total distance from each client
// to its nearest open facility.
//
// Input JSON:
//   { "num_clients": N, "num_facilities": M,
//     "distances": [[d_ij for j in 0..M-1] for i in 0..N-1], "k": K }
//
// Output JSON:
//   { "algorithm": "kmedian", "k": K,
//     "opened_facilities_greedy": [...], "total_cost_greedy": C,
//     "opened_facilities_local_search": [...], "total_cost_local_search": C,
//     "assignments_greedy": [...], "assignments_local_search": [...] }
// =============================================================================

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <random>
#include <iomanip>

#include "json.hpp"

using json = nlohmann::json;

namespace aal {

using Graph = std::map<int, std::map<int, double>>;

struct KMedianInstance {
    int num_clients;
    int num_facilities;
    int k;
    std::vector<std::vector<double>> distances;  // distances[i][j] = dist(client i, facility j)
};

std::pair<std::vector<int>, double> compute_assignment(
    const std::vector<std::vector<double>>& distances,
    const std::set<int>& open_facilities,
    int num_clients
) {
    std::vector<int> assignment(num_clients, -1);
    double total_cost = 0.0;

    for (int i = 0; i < num_clients; ++i) {
        double best_dist = std::numeric_limits<double>::infinity();
        int best_facility = -1;
        for (int j : open_facilities) {
            if (distances[i][j] < best_dist) {
                best_dist = distances[i][j];
                best_facility = j;
            }
        }
        assignment[i] = best_facility;
        total_cost += best_dist;
    }

    return {assignment, total_cost};
}

double assignment_cost(
    const std::vector<std::vector<double>>& distances,
    const std::set<int>& open_facilities,
    int num_clients
) {
    double total = 0.0;
    for (int i = 0; i < num_clients; ++i) {
        double best_dist = std::numeric_limits<double>::infinity();
        for (int j : open_facilities) {
            if (distances[i][j] < best_dist) {
                best_dist = distances[i][j];
            }
        }
        total += best_dist;
    }
    return total;
}

std::pair<std::vector<int>, double> kmedian_greedy(
    const std::vector<std::vector<double>>& distances,
    int num_clients,
    int num_facilities,
    int k
) {
    std::set<int> open_facilities;
    std::set<int> remaining_facilities;
    for (int j = 0; j < num_facilities; ++j) {
        remaining_facilities.insert(j);
    }

    double current_cost = 0.0;
    for (int i = 0; i < num_clients; ++i) {
        double best_dist = std::numeric_limits<double>::infinity();
        for (int j = 0; j < num_facilities; ++j) {
            if (distances[i][j] < best_dist) {
                best_dist = distances[i][j];
            }
        }
        current_cost += best_dist;
    }

    for (int step = 0; step < k; ++step) {
        int best_facility = -1;
        double best_new_cost = current_cost;

        for (int j : remaining_facilities) {
            std::set<int> candidate = open_facilities;
            candidate.insert(j);
            double new_cost = assignment_cost(distances, candidate, num_clients);
            if (new_cost < best_new_cost) {
                best_new_cost = new_cost;
                best_facility = j;
            }
        }

        if (best_facility == -1) break;
        open_facilities.insert(best_facility);
        remaining_facilities.erase(best_facility);
        current_cost = best_new_cost;
    }

    auto [assignment, cost] = compute_assignment(distances, open_facilities, num_clients);
    return {assignment, cost};
}

std::pair<std::set<int>, std::vector<int>, double> kmedian_local_search(
    const std::vector<std::vector<double>>& distances,
    int num_clients,
    int num_facilities,
    int k
) {
    std::set<int> open_facilities;
    for (int j = 0; j < k && j < num_facilities; ++j) {
        open_facilities.insert(j);
    }

    auto [assign, cost] = compute_assignment(distances, open_facilities, num_clients);

    bool improved = true;
    while (improved) {
        improved = false;
        std::set<int> closed;
        for (int j = 0; j < num_facilities; ++j) {
            if (!open_facilities.contains(j)) {
                closed.insert(j);
            }
        }

        for (int o : open_facilities) {
            for (int c : closed) {
                std::set<int> candidate = open_facilities;
                candidate.erase(o);
                candidate.insert(c);
                double new_cost = assignment_cost(distances, candidate, num_clients);
                if (new_cost < cost - 1e-12) {
                    open_facilities = candidate;
                    cost = new_cost;
                    improved = true;
                    break;
                }
            }
            if (improved) break;
        }
    }

    auto [assignment, final_cost] = compute_assignment(distances, open_facilities, num_clients);
    return {open_facilities, assignment, final_cost};
}

std::pair<std::vector<int>, double> kmedian_randomized(
    const std::vector<std::vector<double>>& distances,
    int num_clients,
    int num_facilities,
    int k,
    int num_trials,
    unsigned seed
) {
    std::mt19937 rng(seed);
    std::set<int> best_facilities;
    double best_cost = std::numeric_limits<double>::infinity();

    for (int trial = 0; trial < num_trials; ++trial) {
        std::vector<int> indices(num_facilities);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), rng);

        std::set<int> facilities;
        for (int j = 0; j < k && j < num_facilities; ++j) {
            facilities.insert(indices[j]);
        }

        double cost = assignment_cost(distances, facilities, num_clients);
        if (cost < best_cost) {
            best_cost = cost;
            best_facilities = facilities;
        }
    }

    auto [assignment, final_cost] = compute_assignment(distances, best_facilities, num_clients);
    return {assignment, final_cost};
}

} // namespace aal

void solve(std::istream& in, std::ostream& out) {
    json input;
    in >> input;

    int num_clients = input["num_clients"].get<int>();
    int num_facilities = input["num_facilities"].get<int>();
    int k = input["k"].get<int>();
    auto distances = input["distances"].get<std::vector<std::vector<double>>>();

    using namespace aal;

    auto [greedy_assign, greedy_cost] = kmedian_greedy(distances, num_clients, num_facilities, k);

    std::set<int> greedy_facilities;
    for (int i : greedy_assign) {
        if (i != -1) greedy_facilities.insert(i);
    }

    auto [ls_facilities, ls_assign, ls_cost] = kmedian_local_search(
        distances, num_clients, num_facilities, k
    );

    auto [rand_assign, rand_cost] = kmedian_randomized(
        distances, num_clients, num_facilities, k,
        1000, 42
    );

    std::set<int> rand_facilities;
    for (int i : rand_assign) {
        if (i != -1) rand_facilities.insert(i);
    }

    json output;
    output["algorithm"] = "kmedian";
    output["k"] = k;
    output["num_clients"] = num_clients;
    output["num_facilities"] = num_facilities;

    output["opened_facilities_greedy"] = std::vector<int>(greedy_facilities.begin(), greedy_facilities.end());
    output["total_cost_greedy"] = greedy_cost;
    output["assignments_greedy"] = greedy_assign;

    output["opened_facilities_local_search"] = std::vector<int>(ls_facilities.begin(), ls_facilities.end());
    output["total_cost_local_search"] = ls_cost;
    output["assignments_local_search"] = ls_assign;

    output["opened_facilities_randomized"] = std::vector<int>(rand_facilities.begin(), rand_facilities.end());
    output["total_cost_randomized"] = rand_cost;
    output["assignments_randomized"] = rand_assign;

    out << output.dump(2) << std::endl;
}

int main() {
    solve(std::cin, std::cout);
    return 0;
}

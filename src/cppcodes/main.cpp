#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <functional>
#include "chapters.hpp"

int main() {
    std::cout << "======================================================================\n";
    std::cout << "APPROXIMATION ALGORITHMS WITH C++23\n";
    std::cout << "Companion to Vazirani's 'Approximation Algorithms' (Springer 2001)\n";
    std::cout << "======================================================================\n";

    std::vector<std::pair<std::string, std::function<void()>>> demos = {
        {"Chapter 1: Vertex Cover (2-approx)", demo_vertex_cover},
        {"Chapter 2: Set Cover (H_n-approx)", demo_set_cover},
        {"Chapter 3: Steiner Tree & TSP", demo_steiner_tsp},
        {"Chapter 4: Multiway Cut & k-Cut", demo_multiway_kcut},
        {"Chapter 5: k-Center (2-approx)", demo_kcenter},
        {"Chapter 6: Feedback Vertex Set (2-approx)", demo_feedback_vertex_set},
        {"Chapter 7: Shortest Superstring (3-approx)", demo_shortest_superstring},
        {"Chapter 8: Knapsack FPTAS", demo_knapsack},
        {"Chapter 9: Bin Packing (APTAS)", demo_bin_packing},
        {"Chapter 10: Minimum Makespan Scheduling (PTAS)", demo_makespan},
        {"Chapter 11: Euclidean TSP (PTAS Heuristic)", demo_euclidean_tsp},
        {"Chapters 12-14: LP-Based Algorithms", demo_lp_algorithms},
        {"Chapter 13: Knapsack FPTAS (Chapter 13 Edition)", demo_knapsack_ch13},
        {"Chapter 15: Weighted Vertex Cover via Primal-Dual", demo_weighted_vertex_cover_pd},
        {"Chapter 16: Randomized Rounding (Max-SAT)", demo_randomized_rounding},
        {"Chapter 17: Chernoff Bounds (Set Cover)", demo_chernoff_bounds},
        {"Chapter 18: Semidefinite Programming (Max-Cut)", demo_sdp_max_cut},
        {"Chapter 19: Multiway Cut via LP Rounding", demo_multiway_cut_lp},
        {"Chapter 21: Steiner Forest (2-approx)", demo_steiner_forest},
        {"Chapter 22: Steiner Network (Jain's Iterative Rounding)", demo_steiner_network},
        {"Chapter 23: Feedback Vertex Set via Primal-Dual", demo_primal_dual_fvs},
        {"Chapter 24: Facility Location", demo_facility_location},
        {"Chapter 26: Semidefinite Programming (Max 2-SAT)", demo_max_2sat},
        {"Chapter 28: Counting Problems", demo_counting_problems},
        {"Chapter 30: Multicut in Trees (2-approx)", demo_tree_multicut}
    };

    for (const auto& [name, demo] : demos) {
        std::cout << "\n======================================================================\n";
        std::cout << name << "\n";
        std::cout << "======================================================================\n";
        try {
            demo();
        } catch (const std::exception& e) {
            std::cerr << "Error in " << name << ": " << e.what() << "\n";
        } catch (...) {
            std::cerr << "Unknown error in " << name << "\n";
        }
    }

    std::cout << "\n======================================================================\n";
    std::cout << "ALL DEMOS COMPLETED\n";
    std::cout << "======================================================================\n";
    return 0;
}

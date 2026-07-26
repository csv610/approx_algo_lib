#pragma once

#include <vector>
#include <set>
#include <string>
#include <utility>
#include <map>
#include <memory>
#include <tuple>

namespace aal {

// Graph and Edge types used in standard combinatorics
using Edge = std::pair<int, int>;
using Graph = std::vector<std::vector<int>>;
using WeightedGraph = std::map<int, std::map<int, double>>;

// =====================================================================
// Chapter 1: Vertex Cover
// =====================================================================
std::vector<Edge> maximal_matching(const Graph& graph);
std::set<int> vertex_cover_approx_2(const Graph& graph);
std::set<int> vertex_cover_exact_bruteforce(const Graph& graph);
Graph tight_example_k_n_n(int n);

// =====================================================================
// Chapter 2: Set Cover
// =====================================================================
std::pair<std::vector<int>, double> greedy_set_cover(
    const std::set<int>& universe, 
    const std::map<int, std::set<int>>& sets, 
    const std::map<int, double>& costs
);

// =====================================================================
// Chapter 3: Steiner Tree & TSP
// =====================================================================
std::pair<std::vector<Edge>, double> steiner_tree_2approx(const WeightedGraph& graph, const std::set<int>& terminals);
std::pair<std::vector<int>, double> tsp_2approx_mst(const WeightedGraph& graph);
std::pair<std::vector<int>, double> tsp_christofides_1_5_approx(const WeightedGraph& graph);

// =====================================================================
// Chapter 4: Multiway Cut & k-Cut
// =====================================================================
std::pair<std::set<Edge>, double> multiway_cut_2_2k(const WeightedGraph& graph, const std::set<int>& terminals);
std::pair<std::set<Edge>, double> min_k_cut_2_2k(const WeightedGraph& graph, int k);

// =====================================================================
// Chapter 5: k-Center
// =====================================================================
std::pair<std::set<int>, double> kcenter_parametric_pruning(const WeightedGraph& graph, int k);
std::pair<std::set<int>, double> weighted_kcenter_3approx(const WeightedGraph& graph, const std::map<int, double>& weights, double W);

// =====================================================================
// Chapter 6: Feedback Vertex Set
// =====================================================================
using Multigraph = std::map<int, std::vector<int>>;
std::set<int> feedback_vertex_set_approx(Multigraph graph, std::map<int, double> weights);

// =====================================================================
// Chapter 7: Shortest Superstring
// =====================================================================
std::string greedy_superstring(const std::vector<std::string>& strings);
std::string shortest_superstring_4approx(const std::vector<std::string>& strings);
std::string shortest_superstring_3approx(const std::vector<std::string>& strings);

// =====================================================================
// Chapter 8: Knapsack FPTAS
// =====================================================================
std::pair<std::vector<int>, int> knapsack_dp(const std::vector<int>& weights, const std::vector<int>& values, int capacity);
std::pair<std::vector<int>, int> knapsack_fptas(const std::vector<int>& weights, const std::vector<int>& values, int capacity, double epsilon);

// =====================================================================
// Chapter 9: Bin Packing
// =====================================================================
std::vector<std::vector<double>> bin_packing_aptas(std::vector<double> items, double eps = 0.3, double capacity = 1.0);

// =====================================================================
// Chapter 10: Minimum Makespan Scheduling
// =====================================================================
std::vector<std::vector<double>> makespan_ptas(const std::vector<double>& jobs, int m, double eps = 0.25);

// =====================================================================
// Chapter 11: Euclidean TSP
// =====================================================================
std::pair<std::vector<int>, double> quadtree_tsp(const std::vector<std::pair<double, double>>& points);

// =====================================================================
// Chapters 12-14: LP-Based Algorithms
// =====================================================================
std::pair<std::vector<int>, double> set_cover_lp_rounding(const std::set<int>& universe, const std::map<int, std::set<int>>& sets, const std::map<int, double>& costs);
std::pair<std::vector<int>, double> set_cover_primal_dual(const std::set<int>& universe, const std::map<int, std::set<int>>& sets, const std::map<int, double>& costs);
std::pair<std::set<int>, double> vertex_cover_lp_rounding(const std::map<int, std::map<int, double>>& graph);

// =====================================================================
// Chapter 15: Weighted Vertex Cover via Primal-Dual
// =====================================================================
std::pair<std::set<int>, std::map<std::pair<int, int>, double>> vertex_cover_primal_dual(const std::vector<int>& vertices, const std::vector<std::pair<int, int>>& edges, const std::map<int, double>& weights);

// =====================================================================
// Chapter 16: Randomized Rounding (Max-SAT)
// =====================================================================
std::tuple<std::vector<bool>, double, double> randomized_rounding_max_sat(int n_vars, const std::vector<std::pair<std::set<int>, std::set<int>>>& clauses, const std::vector<double>& weights, const std::vector<double>& y_lp, int trials = 500);

// =====================================================================
// Chapter 17: Chernoff Bounds (Set Cover)
// =====================================================================
std::tuple<std::set<int>, double, bool> set_cover_randomized_rounding(const std::set<int>& universe, const std::map<int, std::set<int>>& sets, const std::map<int, double>& costs, const std::vector<double>& x_lp, double c);

// =====================================================================
// Chapter 18: Semidefinite Programming (Max-Cut)
// =====================================================================
std::vector<std::vector<double>> optimize_max_cut_vectors(int n, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& weights, int dim = 2, double lr = 0.1, int steps = 300);
std::tuple<std::vector<int>, double, double> goemans_williamson_max_cut(int n, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& weights, const std::vector<std::vector<double>>& vectors, int trials = 500);

// =====================================================================
// Chapter 19: Multiway Cut via LP Rounding
// =====================================================================
std::tuple<std::vector<std::vector<double>>, double> solve_multiway_cut_lp(int n, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& costs, const std::vector<int>& terminals);
std::tuple<std::vector<std::pair<int, int>>, double> calinescu_karloff_rabani_rounding(int n, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& costs, const std::vector<int>& terminals, const std::vector<std::vector<double>>& d, int trials = 500);

// =====================================================================
// Chapter 21: Steiner Forest
// =====================================================================
std::vector<std::pair<int, int>> steiner_forest_primal_dual(const std::vector<int>& vertices, const std::vector<std::tuple<int, int, double>>& edges, const std::vector<std::pair<int, int>>& pairs);

// =====================================================================
// Chapter 22: Steiner Network
// =====================================================================
std::vector<std::pair<int, int>> jain_iterative_rounding(int n, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& costs, const std::map<std::pair<int, int>, int>& r);

// =====================================================================
// Chapter 23: Feedback Vertex Set via Primal-Dual
// =====================================================================
std::vector<int> primal_dual_fvs(const std::vector<int>& vertices, const std::vector<std::pair<int, int>>& edges, const std::map<int, double>& weights);

// =====================================================================
// Chapter 24: Facility Location
// =====================================================================
struct Facility { double cost; std::map<int, double> clients; };
std::tuple<std::set<int>, std::map<int, int>, double> facility_location_greedy(const std::map<int, Facility>& facilities, const std::vector<int>& clients);
std::tuple<std::set<int>, std::map<int, int>, double> facility_location_primal_dual(const std::map<int, Facility>& facilities, const std::vector<int>& clients);

// =====================================================================
// Chapter 26: Semidefinite Programming (Max 2-SAT)
// =====================================================================
std::vector<std::vector<double>> optimize_max_2sat_vectors(int n_vars, const std::vector<std::pair<int, int>>& clauses, const std::vector<double>& weights, int dim = 3, double lr = 0.1, int steps = 300);
std::tuple<std::vector<bool>, double, double> goemans_williamson_max_2sat(int n_vars, const std::vector<std::pair<int, int>>& clauses, const std::vector<double>& weights, const std::vector<std::vector<double>>& vectors, int trials = 500);

// =====================================================================
// Chapter 30: Multicut in Trees
// =====================================================================
std::vector<std::pair<int, int>> multicut_in_trees(int n, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& costs, const std::vector<std::pair<int, int>>& pairs);

} // namespace aal

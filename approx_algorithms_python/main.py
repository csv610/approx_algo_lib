"""
Approximation Algorithms with Python - Main Runner
===================================================
Runs all chapters demonstrating algorithms from Vazirani's book.
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from ch01_intro import demo_vertex_cover
from ch02_set_cover import demo_set_cover
from ch03_steiner_tsp import demo_steiner_tsp
from ch04_multiway_kcut import demo_multiway_kcut
from ch05_kcenter import demo_kcenter
from ch06_feedback_vertex_set import demo_feedback_vertex_set
from ch07_shortest_superstring import demo_shortest_superstring
from ch08_knapsack import demo_knapsack
from ch09_bin_packing import demo_bin_packing
from ch10_makespan import demo_makespan
from ch11_euclidean_tsp import demo_euclidean_tsp
from ch12_14_lp_algorithms import demo_lp_algorithms
from ch15_weighted_vertex_cover_pd import demo_weighted_vertex_cover_pd
from ch16_randomized_rounding import demo_randomized_rounding
from ch17_chernoff_bounds import demo_chernoff_bounds
from ch18_sdp_maxcut import demo_sdp_max_cut
from ch19_multiway_rounding import demo_multiway_cut_lp
from ch21_steiner_forest import demo_steiner_forest
from ch22_steiner_network import demo_steiner_network
from ch23_primal_dual_fvs import demo_primal_dual_fvs
from ch24_facility_location import demo_facility_location
from ch26_sdp_sat import demo_max_2sat
from ch30_tree_multicut import demo_tree_multicut


def main():
    print("=" * 70)
    print("APPROXIMATION ALGORITHMS WITH PYTHON")
    print("Companion to Vazirani's 'Approximation Algorithms' (Springer 2001)")
    print("=" * 70)
    
    demos = [
        ("Chapter 1: Vertex Cover (2-approx)", demo_vertex_cover),
        ("Chapter 2: Set Cover (H_n-approx)", demo_set_cover),
        ("Chapter 3: Steiner Tree & TSP", demo_steiner_tsp),
        ("Chapter 4: Multiway Cut & k-Cut", demo_multiway_kcut),
        ("Chapter 5: k-Center (2-approx)", demo_kcenter),
        ("Chapter 6: Feedback Vertex Set (2-approx)", demo_feedback_vertex_set),
        ("Chapter 7: Shortest Superstring (3-approx)", demo_shortest_superstring),
        ("Chapter 8: Knapsack FPTAS", demo_knapsack),
        ("Chapter 9: Bin Packing (APTAS)", demo_bin_packing),
        ("Chapter 10: Minimum Makespan Scheduling (PTAS)", demo_makespan),
        ("Chapter 11: Euclidean TSP (PTAS Heuristic)", demo_euclidean_tsp),
        ("Chapters 12-14: LP-Based Algorithms", demo_lp_algorithms),
        ("Chapter 15: Weighted Vertex Cover via Primal-Dual", demo_weighted_vertex_cover_pd),
        ("Chapter 16: Randomized Rounding (Max-SAT)", demo_randomized_rounding),
        ("Chapter 17: Chernoff Bounds (Set Cover)", demo_chernoff_bounds),
        ("Chapter 18: Semidefinite Programming (Max-Cut)", demo_sdp_max_cut),
        ("Chapter 19: Multiway Cut via LP Rounding", demo_multiway_cut_lp),
        ("Chapter 21: Steiner Forest (2-approx)", demo_steiner_forest),
        ("Chapter 22: Steiner Network (Jain's Iterative Rounding)", demo_steiner_network),
        ("Chapter 23: Feedback Vertex Set via Primal-Dual", demo_primal_dual_fvs),
        ("Chapter 24: Facility Location", demo_facility_location),
        ("Chapter 26: Semidefinite Programming (Max 2-SAT)", demo_max_2sat),
        ("Chapter 30: Multicut in Trees (2-approx)", demo_tree_multicut),
    ]
    
    for name, demo in demos:
        print(f"\n{'=' * 70}")
        print(name)
        print(f"{'=' * 70}")
        try:
            demo()
        except Exception as e:
            print(f"Error in {name}: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 70)
    print("ALL DEMOS COMPLETED")
    print("=" * 70)


if __name__ == "__main__":
    main()
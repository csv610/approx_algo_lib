"""
Approximation Algorithms with Python
=====================================
A concise implementation-focused companion to Vazirani's 
"Approximation Algorithms" (Springer 2001).

Pure Python, zero dependencies.
"""

__version__ = "1.0.0"
__author__ = "Companion to Vazirani (2001)"

# Chapter 1: Introduction (Vertex Cover)
from .intro import (
    vertex_cover_approx_2,
    vertex_cover_approx_2_edge_weighted,
    maximal_matching,
    tight_example_k_n_n,
    vertex_cover_exact_bruteforce,
)

# Chapter 2: Set Cover
from .set_cover import (
    greedy_set_cover,
    set_cover_tight_example,
    harmonic_number,
)

# Chapter 3: Steiner Tree & TSP
from .steiner_tsp import (
    mst_prim,
    steiner_tree_2approx,
    tsp_2approx_mst,
    tsp_christofides_1_5_approx,
    metric_closure,
)

# Chapter 4: Multiway Cut & k-Cut
from .multiway_kcut import (
    multiway_cut_2_2k,
    min_k_cut_2_2k,
    gomory_hu_tree,
)

# Chapter 5: k-Center
from .kcenter import (
    kcenter_parametric_pruning,
    weighted_kcenter_3approx,
)

# Chapter 6: Feedback Vertex Set
from .feedback_vertex_set import (
    feedback_vertex_set_approx,
    is_feedback_vertex_set,
)

# Chapter 7: Shortest Superstring
from .shortest_superstring import (
    shortest_superstring_3approx,
    shortest_superstring_4approx,
    greedy_superstring,
)

# Chapter 8: Knapsack FPTAS
from .knapsack import (
    knapsack_dp,
    knapsack_fptas,
    knapsack_greedy_ratio,
    knapsack_unbounded_dp,
)

# Chapter 9: Bin Packing
from .bin_packing import (
    bin_packing_aptas,
)

# Chapter 10: Minimum Makespan Scheduling
from .makespan import (
    makespan_ptas,
    list_scheduling,
    lpt_scheduling,
)

# Chapter 11: Euclidean TSP
from .euclidean_tsp import (
    quadtree_tsp,
    held_karp_tsp,
)

# Chapters 12-14: LP-Based Algorithms (Set Cover & Vertex Cover)
from .lp_algorithms import (
    set_cover_lp_rounding,
    set_cover_primal_dual,
    vertex_cover_lp_rounding,
)

# Chapter 15: Weighted Vertex Cover via Primal-Dual
from .weighted_vertex_cover_pd import (
    vertex_cover_primal_dual,
)

# Chapter 16: Randomized Rounding (Max-SAT)
from .randomized_rounding import (
    solve_max_sat_lp,
    randomized_rounding_max_sat,
)

# Chapter 17: Chernoff Bounds (Set Cover)
from .chernoff_bounds import (
    solve_set_cover_lp,
    set_cover_randomized_rounding,
)

# Chapter 18: Semidefinite Programming (Max-Cut)
from .sdp_maxcut import (
    goemans_williamson_max_cut,
)

# Chapter 19: Multiway Cut via LP Rounding
from .multiway_rounding import (
    solve_multiway_cut_lp,
    calinescu_karloff_rabani_rounding,
)

# Chapter 21: Steiner Forest
from .steiner_forest import (
    steiner_forest_primal_dual,
)

# Chapter 22: Steiner Network
from .steiner_network import (
    jain_iterative_rounding,
)

# Chapter 23: Feedback Vertex Set via Primal-Dual
from .primal_dual_fvs import (
    primal_dual_fvs,
)

# Chapter 24: Facility Location
from .facility_location import (
    facility_location_greedy,
    facility_location_primal_dual,
    facility_location_lp_rounding,
)

# Chapter 26: Semidefinite Programming (Max 2-SAT)
from .sdp_sat import (
    goemans_williamson_max_2sat,
)

# Chapter 30: Multicut in Trees
from .tree_multicut import (
    multicut_in_trees,
)

__all__ = [
    # Ch 1
    'vertex_cover_approx_2', 'vertex_cover_approx_2_edge_weighted', 'maximal_matching', 
    'tight_example_k_n_n', 'vertex_cover_exact_bruteforce',
    # Ch 2
    'greedy_set_cover', 'set_cover_tight_example', 'harmonic_number',
    # Ch 3
    'mst_prim', 'steiner_tree_2approx', 'tsp_2approx_mst',
    'tsp_christofides_1_5_approx', 'metric_closure',
    # Ch 4
    'multiway_cut_2_2k', 'min_k_cut_2_2k', 'gomory_hu_tree',
    # Ch 5
    'kcenter_parametric_pruning', 'weighted_kcenter_3approx',
    # Ch 6
    'feedback_vertex_set_approx', 'is_feedback_vertex_set',
    # Ch 7
    'shortest_superstring_3approx', 'shortest_superstring_4approx', 'greedy_superstring',
    # Ch 8
    'knapsack_dp', 'knapsack_fptas', 'knapsack_greedy_ratio', 'knapsack_unbounded_dp',
    # Ch 9
    'bin_packing_aptas',
    # Ch 10
    'makespan_ptas', 'list_scheduling', 'lpt_scheduling',
    # Ch 11
    'quadtree_tsp', 'held_karp_tsp',
    # Ch 12-14
    'set_cover_lp_rounding', 'set_cover_primal_dual', 'vertex_cover_lp_rounding',
    # Ch 15
    'vertex_cover_primal_dual',
    # Ch 16
    'solve_max_sat_lp', 'randomized_rounding_max_sat',
    # Ch 17
    'solve_set_cover_lp', 'set_cover_randomized_rounding',
    # Ch 18
    'goemans_williamson_max_cut',
    # Ch 19
    'solve_multiway_cut_lp', 'calinescu_karloff_rabani_rounding',
    # Ch 21
    'steiner_forest_primal_dual',
    # Ch 22
    'jain_iterative_rounding',
    # Ch 23
    'primal_dual_fvs',
    # Ch 24
    'facility_location_greedy', 'facility_location_primal_dual', 'facility_location_lp_rounding',
    # Ch 26
    'goemans_williamson_max_2sat',
    # Ch 30
    'multicut_in_trees',
]
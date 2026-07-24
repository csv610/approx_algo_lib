# Approximation Algorithms with Python

A concise, implementation-focused companion to **Vijay Vazirani's "Approximation Algorithms" (Springer 2001)**.

This project provides clean, readable Python implementations of the core approximation algorithms from the book, with zero external dependencies.

## Structure

| Chapter | Algorithm | Approximation | File |
|---------|-----------|---------------|------|
| 1 | Vertex Cover (maximal matching) | 2 | `ch01_intro.py` |
| 2 | Set Cover (greedy) | Hₙ | `ch02_set_cover.py` |
| 3 | Metric Steiner Tree / TSP | 2, 3/2 | `ch03_steiner_tsp.py` |
| 4 | Multiway Cut / k-Cut | 2-2/k | `ch04_multiway_kcut.py` |
| 5 | k-Center (parametric pruning) | 2, 3 | `ch05_kcenter.py` |
| 6 | Feedback Vertex Set | 2 | `ch06_feedback_vertex_set.py` |
| 7 | Shortest Superstring | 3 | `ch07_shortest_superstring.py` |
| 8 | Knapsack FPTAS | (1-ε) | `ch08_knapsack.py` |
| 9 | Bin Packing (APTAS) | (1+ε) | `ch09_bin_packing.py` |
| 10 | Minimum Makespan Scheduling | (1+ε) | `ch10_makespan.py` |
| 11 | Euclidean TSP (PTAS Heuristic) | PTAS | `ch11_euclidean_tsp.py` |
| 12-14 | LP Rounding / Primal-Dual | f | `ch12_14_lp_algorithms.py` |
| 13 | Knapsack FPTAS (Duplicate) | (1-ε) | `ch13_knapsack.py` |
| 15 | Weighted Vertex Cover via Primal-Dual | 2 | `ch15_weighted_vertex_cover_pd.py` |
| 16 | Randomized Rounding (Max-SAT) | (1-1/e) | `ch16_randomized_rounding.py` |
| 17 | Chernoff Bounds (Set Cover) | O(log n) | `ch17_chernoff_bounds.py` |
| 18 | Semidefinite Programming (Max-Cut) | 0.878 | `ch18_sdp_maxcut.py` |
| 19 | Multiway Cut LP Rounding (CKR) | 1.5 | `ch19_multiway_rounding.py` |
| 21 | Steiner Forest | 2 | `ch21_steiner_forest.py` |
| 22 | Steiner Network | 2 | `ch22_steiner_network.py` |
| 23 | Feedback Vertex Set via Primal-Dual | 2 | `ch23_primal_dual_fvs.py` |
| 24 | Facility Location | 3 | `ch24_facility_location.py` |
| 26 | Max 2-SAT (SDP) | 0.878 | `ch26_sdp_sat.py` |
| 30 | Multicut in Trees | 2 | `ch30_tree_multicut.py` |

## Running

```bash
cd approx_algorithms_python
python3 main.py
```

## Example Output

```
======================================================================
Chapter 1: Vertex Cover - Factor 2 Approximation
======================================================================

1. Tight Example: K_{4,4}
  Graph: K_{4,4} (8 vertices, 16 edges)
  Approx cover size: 8
  Optimal cover size: 4
  Ratio: 2.00
  (Optimal picks one side: 4 vertices)

2. Random Graph (10 vertices, p=0.3)
  Graph: 10 vertices, 17 edges
  Approx cover size: 8
  Optimal cover size: 6
  Ratio: 1.33
```

## Key Features

- **Pure Python** - no dependencies (only stdlib)
- **Educational** - each algorithm includes the Vazirani chapter reference
- **Runnable demos** - every chapter has `demo_*()` showing tight examples
- **Zero-setup** - runs immediately with standard Python 3.8+

## References

- Vazirani, V. V. *Approximation Algorithms*. Springer, 2001.
- PAAL Library: https://paal.mimuw.edu.pl (C++ reference implementations)

## License

MIT - Free for educational and commercial use.
# Approximation Algorithms with Python - Summary

## Project Structure
```
/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/
├── __init__.py                    # Package exports
├── main.py                        # Runs all demos
├── README.md                      # Documentation
├── book.py                        # Package info
├── ch01_intro.py                  # Ch 1: Vertex Cover (2-approx)
├── ch02_set_cover.py              # Ch 2: Set Cover (H_n-approx)
├── ch03_steiner_tsp.py            # Ch 3: Steiner Tree & TSP (2, 3/2-approx)
├── ch04_multiway_kcut.py          # Ch 4: Multiway Cut & k-Cut (2-2/k)
├── ch05_kcenter.py                # Ch 5: k-Center (2, 3-approx)
├── ch06_feedback_vertex_set.py    # Ch 6: Feedback Vertex Set (2-approx)
├── ch07_shortest_superstring.py    # Ch 7: Shortest Superstring (3-approx)
├── ch08_knapsack.py               # Ch 8: Knapsack FPTAS
├── ch09_bin_packing.py            # Ch 9: Bin Packing (APTAS)
├── ch10_makespan.py               # Ch 10: Minimum Makespan Scheduling (PTAS)
├── ch11_euclidean_tsp.py          # Ch 11: Euclidean TSP (PTAS Heuristic)
├── ch12_14_lp_algorithms.py       # Ch 12-14: LP Rounding & Primal-Dual
├── ch13_knapsack.py               # Ch 13: Knapsack FPTAS (Duplicate of Ch 8)
├── ch15_weighted_vertex_cover_pd.py # Ch 15: Weighted Vertex Cover via Primal-Dual (2-approx)
├── ch16_randomized_rounding.py    # Ch 16: Randomized Rounding (Max-SAT) (1-1/e-approx)
├── ch17_chernoff_bounds.py        # Ch 17: Chernoff Bounds (Set Cover)
├── ch18_sdp_maxcut.py             # Ch 18: Semidefinite Programming (Max-Cut) (0.878-approx)
├── ch19_multiway_rounding.py      # Ch 19: Multiway Cut via LP Rounding (1.5-approx)
├── ch21_steiner_forest.py         # Ch 21: Steiner Forest (2-approx)
├── ch22_steiner_network.py        # Ch 22: Steiner Network (Jain's Iterative Rounding) (2-approx)
├── ch23_primal_dual_fvs.py        # Ch 23: Feedback Vertex Set via Primal-Dual (2-approx)
├── ch24_facility_location.py      # Ch 24: Facility Location (3-approx)
├── ch26_sdp_sat.py                # Ch 26: Semidefinite Programming (Max 2-SAT) (0.878-approx)
└── ch30_tree_multicut.py          # Ch 30: Multicut in Trees (2-approx)
```

## Algorithms Implemented

### Part I: Combinatorial Algorithms (Ch 1-11)
| Algorithm | Technique | Approximation | Key File |
|-----------|-----------|---------------|----------|
| Vertex Cover | Maximal matching | 2 | ch01_intro.py |
| Set Cover | Greedy | Hₙ | ch02_set_cover.py |
| Metric Steiner Tree | Metric closure + MST | 2 | ch03_steiner_tsp.py |
| Metric TSP | MST double-tree | 2 | ch03_steiner_tsp.py |
| Metric TSP | Christofides | 3/2 | ch03_steiner_tsp.py |
| Multiway Cut | Isolating cuts | 2-2/k | ch04_multiway_kcut.py |
| Minimum k-Cut | Gomory-Hu tree | 2-2/k | ch04_multiway_kcut.py |
| k-Center | Parametric pruning | 2 | ch05_kcenter.py |
| Weighted k-Center | Parametric pruning | 3 | ch05_kcenter.py |
| Feedback Vertex Set | Layering / Local ratio | 2 | ch06_feedback_vertex_set.py |
| Shortest Superstring | Cycle cover greedy | 3 | ch07_shortest_superstring.py |
| Knapsack | DP with rounding | (1-ε) | ch08_knapsack.py |
| Bin Packing | APTAS | (1+ε) | ch09_bin_packing.py |
| Minimum Makespan | PTAS | (1+ε) | ch10_makespan.py |
| Euclidean TSP | Quadtree heuristic | PTAS | ch11_euclidean_tsp.py |

### Part II: LP & SDP-Based Algorithms (Ch 12-30)
| Algorithm | Technique | Approximation | Key File |
|-----------|-----------|---------------|----------|
| Set Cover | LP Rounding | f | ch12_14_lp_algorithms.py |
| Set Cover | Primal-Dual | f | ch12_14_lp_algorithms.py |
| Vertex Cover | LP Rounding | 2 | ch12_14_lp_algorithms.py |
| Knapsack | FPTAS (value scaling) | (1-ε) | ch13_knapsack.py |
| Weighted Vertex Cover | Primal-Dual | 2 | ch15_weighted_vertex_cover_pd.py |
| Max-SAT | Randomized rounding | (1-1/e) | ch16_randomized_rounding.py |
| Set Cover | Chernoff randomized rounding | O(log n) | ch17_chernoff_bounds.py |
| Max-Cut | SDP + Hyperplane rounding | 0.878 | ch18_sdp_maxcut.py |
| Multiway Cut | LP Rounding (CKR) | 1.5 | ch19_multiway_rounding.py |
| Steiner Forest | Primal-Dual (AKR) | 2 | ch21_steiner_forest.py |
| Steiner Network | Jain's Iterative rounding | 2 | ch22_steiner_network.py |
| Feedback Vertex Set | Primal-Dual | 2 | ch23_primal_dual_fvs.py |
| Facility Location | Greedy / Primal-Dual / Local Search | 3 | ch24_facility_location.py |
| Max 2-SAT | SDP + Hyperplane rounding | 0.878 | ch26_sdp_sat.py |
| Multicut in Trees | Primal-Dual | 2 | ch30_tree_multicut.py |

## Running
```bash
cd /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python
python3 main.py
```

## Key Features
- **Pure Python** - zero dependencies (stdlib only)
- **Educational** - each algorithm references Vazirani chapter
- **Runnable demos** - every chapter includes tight examples
- **Clean implementations** - focus on readability over optimization
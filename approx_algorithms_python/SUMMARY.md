# Approximation Algorithms with Python - Summary

## Project Structure
```
/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/
├── __init__.py                    # Package exports
├── main.py                        # Runs all demos
├── README.md                      # Documentation
├── book.py                        # Package info
├── intro.py                  # Ch 1: Vertex Cover (2-approx)
├── set_cover.py              # Ch 2: Set Cover (H_n-approx)
├── steiner_tsp.py            # Ch 3: Steiner Tree & TSP (2, 3/2-approx)
├── multiway_kcut.py          # Ch 4: Multiway Cut & k-Cut (2-2/k)
├── kcenter.py                # Ch 5: k-Center (2, 3-approx)
├── feedback_vertex_set.py    # Ch 6: Feedback Vertex Set (2-approx)
├── shortest_superstring.py    # Ch 7: Shortest Superstring (3-approx)
├── knapsack.py               # Ch 8: Knapsack FPTAS
├── bin_packing.py            # Ch 9: Bin Packing (APTAS)
├── makespan.py               # Ch 10: Minimum Makespan Scheduling (PTAS)
├── euclidean_tsp.py          # Ch 11: Euclidean TSP (PTAS Heuristic)
├── lp_algorithms.py       # Ch 12-14: LP Rounding & Primal-Dual
├── knapsack_ch13.py               # Ch 13: Knapsack FPTAS (Duplicate of Ch 8)
├── weighted_vertex_cover_pd.py # Ch 15: Weighted Vertex Cover via Primal-Dual (2-approx)
├── randomized_rounding.py    # Ch 16: Randomized Rounding (Max-SAT) (1-1/e-approx)
├── chernoff_bounds.py        # Ch 17: Chernoff Bounds (Set Cover)
├── sdp_maxcut.py             # Ch 18: Semidefinite Programming (Max-Cut) (0.878-approx)
├── multiway_rounding.py      # Ch 19: Multiway Cut via LP Rounding (1.5-approx)
├── steiner_forest.py         # Ch 21: Steiner Forest (2-approx)
├── steiner_network.py        # Ch 22: Steiner Network (Jain's Iterative Rounding) (2-approx)
├── primal_dual_fvs.py        # Ch 23: Feedback Vertex Set via Primal-Dual (2-approx)
├── facility_location.py      # Ch 24: Facility Location (3-approx)
├── sdp_sat.py                # Ch 26: Semidefinite Programming (Max 2-SAT) (0.878-approx)
└── tree_multicut.py          # Ch 30: Multicut in Trees (2-approx)
```

## Algorithms Implemented

### Part I: Combinatorial Algorithms (Ch 1-11)
| Algorithm | Technique | Approximation | Key File |
|-----------|-----------|---------------|----------|
| Vertex Cover | Maximal matching | 2 | intro.py |
| Set Cover | Greedy | Hₙ | set_cover.py |
| Metric Steiner Tree | Metric closure + MST | 2 | steiner_tsp.py |
| Metric TSP | MST double-tree | 2 | steiner_tsp.py |
| Metric TSP | Christofides | 3/2 | steiner_tsp.py |
| Multiway Cut | Isolating cuts | 2-2/k | multiway_kcut.py |
| Minimum k-Cut | Gomory-Hu tree | 2-2/k | multiway_kcut.py |
| k-Center | Parametric pruning | 2 | kcenter.py |
| Weighted k-Center | Parametric pruning | 3 | kcenter.py |
| Feedback Vertex Set | Layering / Local ratio | 2 | feedback_vertex_set.py |
| Shortest Superstring | Cycle cover greedy | 3 | shortest_superstring.py |
| Knapsack | DP with rounding | (1-ε) | knapsack.py |
| Bin Packing | APTAS | (1+ε) | bin_packing.py |
| Minimum Makespan | PTAS | (1+ε) | makespan.py |
| Euclidean TSP | Quadtree heuristic | PTAS | euclidean_tsp.py |

### Part II: LP & SDP-Based Algorithms (Ch 12-30)
| Algorithm | Technique | Approximation | Key File |
|-----------|-----------|---------------|----------|
| Set Cover | LP Rounding | f | lp_algorithms.py |
| Set Cover | Primal-Dual | f | lp_algorithms.py |
| Vertex Cover | LP Rounding | 2 | lp_algorithms.py |
| Knapsack | FPTAS (value scaling) | (1-ε) | knapsack_ch13.py |
| Weighted Vertex Cover | Primal-Dual | 2 | weighted_vertex_cover_pd.py |
| Max-SAT | Randomized rounding | (1-1/e) | randomized_rounding.py |
| Set Cover | Chernoff randomized rounding | O(log n) | chernoff_bounds.py |
| Max-Cut | SDP + Hyperplane rounding | 0.878 | sdp_maxcut.py |
| Multiway Cut | LP Rounding (CKR) | 1.5 | multiway_rounding.py |
| Steiner Forest | Primal-Dual (AKR) | 2 | steiner_forest.py |
| Steiner Network | Jain's Iterative rounding | 2 | steiner_network.py |
| Feedback Vertex Set | Primal-Dual | 2 | primal_dual_fvs.py |
| Facility Location | Greedy / Primal-Dual / Local Search | 3 | facility_location.py |
| Max 2-SAT | SDP + Hyperplane rounding | 0.878 | sdp_sat.py |
| Multicut in Trees | Primal-Dual | 2 | tree_multicut.py |

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
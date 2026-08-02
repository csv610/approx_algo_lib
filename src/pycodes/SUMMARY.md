# Approximation Algorithms with Python - Summary

## Project Structure
```
src/pycodes/
├── __init__.py                    # Package exports
├── main.py                        # Runs all demos
├── book.py                        # Package info
├── intro.py                       # Ch 1: Vertex Cover (2-approx)
├── set_cover.py                   # Ch 2: Set Cover (H_n-approx)
├── steiner_tsp.py                 # Ch 3: Steiner Tree & TSP (2, 3/2-approx)
├── multiway_kcut.py               # Ch 4: Multiway Cut & k-Cut (2-2/k)
├── kcenter.py                     # Ch 5: k-Center (2, 3-approx)
├── feedback_vertex_set.py         # Ch 6: Feedback Vertex Set (2-approx)
├── shortest_superstring.py        # Ch 7: Shortest Superstring (3-approx)
├── knapsack.py                    # Ch 8: Knapsack FPTAS
├── bin_packing.py                 # Ch 9: Bin Packing (APTAS)
├── makespan.py                    # Ch 10: Minimum Makespan Scheduling (PTAS)
├── euclidean_tsp.py               # Ch 11: Euclidean TSP (PTAS Heuristic)
├── lp_algorithms.py               # Ch 12-14: LP Rounding & Primal-Dual
├── knapsack_ch13.py               # Ch 13: Knapsack FPTAS (Duplicate)
├── weighted_vertex_cover_pd.py    # Ch 15: Weighted Vertex Cover via PD (2-approx)
├── randomized_rounding.py         # Ch 16: Randomized Rounding (Max-SAT)
├── chernoff_bounds.py             # Ch 17: Chernoff Bounds (Set Cover)
├── sdp_maxcut.py                  # Ch 18: SDP (Max-Cut) (0.878-approx)
├── multiway_rounding.py           # Ch 19: Multiway Cut via LP Rounding (1.5-approx)
├── steiner_forest.py              # Ch 21: Steiner Forest (2-approx)
├── steiner_network.py             # Ch 22: Steiner Network (Jain's) (2-approx)
├── primal_dual_fvs.py             # Ch 23: FVS via Primal-Dual (2-approx)
├── facility_location.py           # Ch 24: Facility Location (3-approx)
├── kmedian.py                     # Ch 25: k-Median (O(log k), 3+eps)
├── sdp_sat.py                     # Ch 26: SDP (Max 2-SAT) (0.878-approx)
├── shortest_vector.py             # Ch 27: Shortest Vector (Gauss, LLL)
├── counting_problems.py           # Ch 28: Counting (DNF, Reliability)
├── hardness_approx.py             # Ch 29: Hardness of Approximation (PCP)
├── tree_multicut.py               # Ch 30: Multicut in Trees (2-approx)
├── unrelated_scheduling.py        # Unrelated Scheduling (2-approx)
├── multicut_general.py            # Multicut General (O(log k)-approx)
├── sparsest_cut.py                # Sparsest Cut (O(sqrt(log n))-approx)
├── multicut_trees.py              # Multicut in Trees (exact)
│
│  --- Williamson & Shmoys Algorithms ---
│
├── edge_coloring.py               # Ch 2.7: Edge Coloring (Delta+1 exact)
├── prize_collecting_steiner.py    # Ch 4.5: Prize-Collecting Steiner (3-approx)
├── shortest_st_path_pd.py         # Ch 7.3: Shortest s-t Path via PD (exact)
├── generalized_steiner.py         # Ch 7.4: Generalized Steiner Tree (2-approx)
├── min_knapsack_pd.py             # Ch 7.5: Min Knapsack via PD (2-approx)
├── balanced_cuts.py               # Ch 8.4: Balanced Cuts (greedy + contraction)
├── tree_metrics.py                # Ch 8.5: FRT Tree Metric Embedding (O(log n))
├── buy_at_bulk.py                 # Ch 8.6: Buy-at-Bulk Network Design
├── facility_location_ls.py        # Ch 9.1: Facility Location Local Search (3-approx)
├── kmedian_ls.py                  # Ch 9.2: k-Median Local Search (3+eps)
├── iterated_rounding.py           # Ch 11.3: Survivable Network Design (2-approx)
├── steiner_tree_rr.py             # Ch 12.3: Steiner Tree Randomized Rounding
├── rent_or_buy.py                 # Ch 12.2: Rent-or-Buy (3-approx)
└── oblivious_routing.py           # Ch 15.2: Oblivious Routing
```

## Algorithms Implemented

### Vazirani - Part I: Combinatorial Algorithms (Ch 1-11)
| Algorithm | Technique | Approximation | File |
|-----------|-----------|---------------|------|
| Vertex Cover | Maximal matching | 2 | intro.py |
| Set Cover | Greedy | H_n | set_cover.py |
| Metric Steiner Tree | Metric closure + MST | 2 | steiner_tsp.py |
| Metric TSP | MST double-tree | 2 | steiner_tsp.py |
| Metric TSP | Christofides | 3/2 | steiner_tsp.py |
| Multiway Cut | Isolating cuts | 2-2/k | multiway_kcut.py |
| Minimum k-Cut | Gomory-Hu tree | 2-2/k | multiway_kcut.py |
| k-Center | Parametric pruning | 2 | kcenter.py |
| Weighted k-Center | Parametric pruning | 3 | kcenter.py |
| Feedback Vertex Set | Layering / Local ratio | 2 | feedback_vertex_set.py |
| Shortest Superstring | Cycle cover greedy | 3 | shortest_superstring.py |
| Knapsack | DP with rounding | (1-e) | knapsack.py |
| Bin Packing | APTAS | (1+e) | bin_packing.py |
| Minimum Makespan | PTAS | (1+e) | makespan.py |
| Euclidean TSP | Quadtree heuristic | PTAS | euclidean_tsp.py |

### Vazirani - Part II: LP & SDP-Based Algorithms (Ch 12-30)
| Algorithm | Technique | Approximation | File |
|-----------|-----------|---------------|------|
| Set Cover | LP Rounding | f | lp_algorithms.py |
| Set Cover | Primal-Dual | f | lp_algorithms.py |
| Vertex Cover | LP Rounding | 2 | lp_algorithms.py |
| Weighted Vertex Cover | Primal-Dual | 2 | weighted_vertex_cover_pd.py |
| Max-SAT | Randomized rounding | (1-1/e) | randomized_rounding.py |
| Set Cover | Chernoff randomized rounding | O(log n) | chernoff_bounds.py |
| Max-Cut | SDP + Hyperplane rounding | 0.878 | sdp_maxcut.py |
| Multiway Cut | LP Rounding (CKR) | 1.5 | multiway_rounding.py |
| Steiner Forest | Primal-Dual (AKR) | 2 | steiner_forest.py |
| Steiner Network | Jain's Iterative rounding | 2 | steiner_network.py |
| Feedback Vertex Set | Primal-Dual | 2 | primal_dual_fvs.py |
| Facility Location | Greedy / PD / Local Search | 3 | facility_location.py |
| Max 2-SAT | SDP + Hyperplane rounding | 0.878 | sdp_sat.py |
| k-Median | LP Rounding + Local Search | O(log k), 3+e | kmedian.py |
| Shortest Vector | Gauss, LLL | 2^{O(n)} | shortest_vector.py |
| DNF Counting | Karp-Luby Monte Carlo | (1+/-e) | counting_problems.py |
| Network Reliability | Monte Carlo | (1+/-e) | counting_problems.py |
| Multicut in Trees | Primal-Dual | 2 | tree_multicut.py |
| Unrelated Scheduling | LP + Parametric Pruning | 2 | unrelated_scheduling.py |
| Multicut (General) | LP + Region Growing | O(log k) | multicut_general.py |
| Sparsest Cut | LP + l1-embedding | O(sqrt(log n)) | sparsest_cut.py |

### Williamson & Shmoys Algorithms
| Algorithm | Chapter | Technique | Approximation | File |
|-----------|---------|-----------|---------------|------|
| Edge Coloring | 2.7 | Greedy + fan recoloring | Delta+1 (exact) | edge_coloring.py |
| Prize-Collecting Steiner | 4.5 | LP rounding | 3 | prize_collecting_steiner.py |
| Shortest s-t Path | 7.3 | Primal-Dual | exact | shortest_st_path_pd.py |
| Generalized Steiner Tree | 7.4 | PD + reverse deletion | 2 | generalized_steiner.py |
| Minimum Knapsack | 7.5 | PD + strengthened LP | 2 | min_knapsack_pd.py |
| Balanced Cuts | 8.4 | Greedy + random contraction | O(1) | balanced_cuts.py |
| FRT Tree Metrics | 8.5 | Recursive clustering | O(log n) dist | tree_metrics.py |
| Buy-at-Bulk | 8.6 | Tree metric embedding | O(log k) | buy_at_bulk.py |
| Facility Location (LS) | 9.1 | Local search | 3 | facility_location_ls.py |
| k-Median (LS) | 9.2 | Local search | 3+e | kmedian_ls.py |
| Survivable Network | 11.3 | Iterated rounding | 2 | iterated_rounding.py |
| Rent-or-Buy | 12.2 | Combined buy/rent | 3 | rent_or_buy.py |
| Steiner Tree (RR) | 12.3 | Randomized rounding | O(log k) | steiner_tree_rr.py |
| Oblivious Routing | 15.2 | Cut-tree packing | competitive | oblivious_routing.py |

## Running
```bash
cd /Users/csv610/Projects/MyBooks/ApproxAlgo/src/pycodes
python3 main.py          # Run all demos
python3 edge_coloring.py # Run specific algorithm
```

## Key Features
- **Pure Python** - zero dependencies (stdlib only)
- **Educational** - each algorithm references Williamson & Shmoys chapter
- **Runnable demos** - every chapter includes tight examples
- **Clean implementations** - focus on readability over optimization

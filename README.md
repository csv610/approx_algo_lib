# Approximation Algorithms: Companion Code for Vazirani's Book

This repository contains complete, runnable implementation companions for Vijay Vazirani’s seminal textbook, **"Approximation Algorithms" (Springer, 2001)**. 

It provides dual codebases (Python & C++23) to serve both educational and software engineering purposes, mapping mathematical proofs directly to executable computer code.

---

## 📂 Repository Structure

The repository is divided into two primary codebases:

* **[approx_algorithms_python/](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/)**:
  - Contains Python implementations of 24 algorithm chapters.
  - Serves as highly readable, interactive "executable pseudocode".
  - Includes the complete LaTeX book file **`book.tex`** and compiled PDF **`book.pdf`** (including the new Chapter 19 LP rounding formulations).
* **[approx_algorithms_cpp/](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/)**:
  - A modern, type-safe, static C++23 library named **`aal`** (Approximate Algorithm Library).
  - Designed with **zero external dependencies** (includes from-scratch Simplex and SDP solvers).
  - Includes a unified public API header **`approx_algorithms.hpp`** and **23 standalone, compile-ready example files** under the `examples/` directory.

---

## 🚀 Quick Start

### 🐍 Python Codebase
To run all 24 Python demos:
```bash
cd approx_algorithms_python
python3 main.py
```

### 🛠️ C++23 aal Library & Examples
To configure, compile, and build the static library, the main runner, and all 23 standalone example binaries:
```bash
cd approx_algorithms_cpp
cmake -B build -S .
cmake --build build
```
Once compiled, you can run the main test runner or any individual example:
```bash
./build/approx_algo_demos
./build/examples/ex_01_vertex_cover
./build/examples/ex_03_steiner_tsp
./build/examples/ex_13_knapsack
```
For integration details, read the C++ **[USAGE.md](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/USAGE.md)**.

---

## 📖 Implemented Chapters Index

This repository covers 24 key algorithmic chapters from Vazirani's book:

| Chapter | Algorithm / Topic | C++ Implementation | Python Implementation |
| :---: | :--- | :---: | :---: |
| **1** | Vertex Cover (2-approx matching) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch01_intro.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch01_intro.py) |
| **2** | Set Cover (Greedy $H_n$-approx) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch02_set_cover.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch02_set_cover.py) |
| **3** | Steiner Tree & TSP (Christofides 1.5-approx) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch03_steiner_tsp.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch03_steiner_tsp.py) |
| **4** | Multiway Cut & $k$-Cut | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch04_multiway_kcut.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch04_multiway_kcut.py) |
| **5** | $k$-Center (Parametric pruning) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch05_kcenter.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch05_kcenter.py) |
| **6** | Feedback Vertex Set (Local Ratio) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch06_feedback_vertex_set.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch06_feedback_vertex_set.py) |
| **7** | Shortest Superstring (Cycle covers) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch07_shortest_superstring.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch07_shortest_superstring.py) |
| **8** | Knapsack (FPTAS) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch08_knapsack.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch08_knapsack.py) |
| **9** | Bin Packing (APTAS) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch09_bin_packing.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch09_bin_packing.py) |
| **10** | Minimum Makespan Scheduling (PTAS) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch10_makespan.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch10_makespan.py) |
| **11** | Euclidean TSP (PTAS Quadtree Partition) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch11_euclidean_tsp.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch11_euclidean_tsp.py) |
| **12-14** | LP algorithms (Rounding, Dual Fitting, Set Cover PD) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch12_14_lp_algorithms.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch12_14_lp_algorithms.py) |
| **13** | Knapsack FPTAS via rounding | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch13_knapsack.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch13_knapsack.py) |
| **15** | Weighted Vertex Cover via Primal-Dual | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch15_weighted_vertex_cover_pd.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch15_weighted_vertex_cover_pd.py) |
| **16** | Randomized Rounding (Max-SAT) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch16_randomized_rounding.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch16_randomized_rounding.py) |
| **17** | Chernoff Bounds (Set Cover rounding) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch17_chernoff_bounds.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch17_chernoff_bounds.py) |
| **18** | SDP for Max-Cut (Goemans-Williamson) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch18_sdp_maxcut.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch18_sdp_maxcut.py) |
| **19** | Multiway Cut LP Rounding (CKR Simplex Rounding) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch19_multiway_rounding.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch19_multiway_rounding.py) |
| **21** | Steiner Forest (PD Schema) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch21_steiner_forest.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch21_steiner_forest.py) |
| **22** | Steiner Network (Jain's Iterative Rounding) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch22_steiner_network.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch22_steiner_network.py) |
| **23** | Feedback Vertex Set via PD | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch23_primal_dual_fvs.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch23_primal_dual_fvs.py) |
| **24** | Facility Location (Greedy, PD, $k$-Median Local Search) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch24_facility_location.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch24_facility_location.py) |
| **26** | SDP for Max 2-SAT | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch26_sdp_sat.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch26_sdp_sat.py) |
| **30** | Multicut in Trees (LCA Depth Pruning) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/ch30_tree_multicut.cpp) | [Source](file:///Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_python/ch30_tree_multicut.py) |

---

## 📜 License
This project is open-source and available under the MIT License.

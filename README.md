# Approximation Algorithms: Companion Code for Vazirani's Book

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/compiler_support/23)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains complete, runnable implementation companions for Vijay Vazirani’s seminal textbook, **"Approximation Algorithms" (Springer, 2001)**. 

It provides dual codebases (Python & C++23) to serve both educational and software engineering purposes, mapping mathematical proofs directly to executable computer code. The C++23 implementations are package-ready and backed by **Google OR-Tools** for industrial-strength Linear Programming.

---

## 🌟 Key Features

* **Dual Implementations**:
  * **Python Codebase**: Highly readable, serving as interactive "executable pseudocode".
  * **C++23 Static Library (`aal`)**: A type-safe, modular implementation under `namespace aal` utilizing modern standard features (`std::print`, ranges, structured bindings).
* **Production-Grade Solvers**:
  - The C++ library integrates **Google OR-Tools (GLOP)** for solving LP relaxations efficiently (replacing standard textbook simplex implementations).
  - Includes a custom spherical gradient-projection Semidefinite Programming (SDP) optimizer built from scratch.
* **Zero Compilation Issues**: Highly portable, zero-warning builds with CMake.
* **100% Chapter Parity**: 24 implemented chapters containing exact algorithms and standalone example binaries.

---

## 📂 Repository Directory Layout

```text
approx_algo_lib/
├── approx_algorithms_python/        # Python Codebase
│   ├── book.tex                    # LaTeX source for the companion notes
│   ├── book.pdf                    # Compiled PDF version of the companion book
│   ├── main.py                     # Main test suite runner for Python
│   └── ch01_intro.py ...           # 24 Chapter Python scripts
├── approx_algorithms_cpp/           # C++23 Codebase
│   ├── CMakeLists.txt              # Build configuration
│   ├── approx_algorithms.hpp       # Unified public API header interface
│   ├── simplex.hpp                 # Google OR-Tools GLOP solver wrapper
│   ├── main.cpp                    # Main test suite runner for C++
│   ├── USAGE.md                    # Integration guide & code reference
│   ├── examples/                   # 23 Standalone, compile-ready example files
│   │   ├── CMakeLists.txt
│   │   ├── ex_01_vertex_cover.cpp
│   │   └── ... ex_30_tree_multicut.cpp
│   └── ch01_intro.cpp ...          # 24 Chapter C++ source implementations
├── README.md                       # Main repository overview
└── .gitignore                      # Residual compile & environment exclusions
```

---

## 🚀 Quick Start & Build Guide

### 🐍 Running Python Demos
Run all 24 Python demos synchronously:
```bash
cd approx_algorithms_python
python3 main.py
```

### 🛠️ Building the C++ Library & Examples

#### Prerequisite: Google OR-Tools
Ensure Google OR-Tools is installed on your system:
* **macOS (Homebrew)**: `brew install or-tools`
* **Ubuntu/Debian**: `sudo apt-get install libortools-dev`

#### Building with CMake
```bash
cd approx_algorithms_cpp
cmake -B build -S .
cmake --build build
```
This compiles and creates:
* **`libapprox_algo.a`** (Static Library target `approx_algo`)
* **`approx_algo_demos`** (Global test runner executable `./build/approx_algo_demos`)
* **23 Standalone Example Binaries** (under `./build/examples/`)

---

## 📖 Implemented Chapters & Approximation Bounds

The repository covers 24 algorithmic chapters from Vazirani's book:

| Chapter | Algorithm / Topic | Approximation Bound | C++ Source | Python Source |
| :---: | :--- | :---: | :---: | :---: |
| **1** | Vertex Cover | Factor-2 (Maximal Matching) | [ch01_intro.cpp](approx_algorithms_cpp/ch01_intro.cpp) | [ch01_intro.py](approx_algorithms_python/ch01_intro.py) |
| **2** | Set Cover | Greedy $H_n$-approx | [ch02_set_cover.cpp](approx_algorithms_cpp/ch02_set_cover.cpp) | [ch02_set_cover.py](approx_algorithms_python/ch02_set_cover.py) |
| **3** | Steiner Tree & TSP | 1.5-approx (Christofides) | [ch03_steiner_tsp.cpp](approx_algorithms_cpp/ch03_steiner_tsp.cpp) | [ch03_steiner_tsp.py](approx_algorithms_python/ch03_steiner_tsp.py) |
| **4** | Multiway Cut & $k$-Cut | $2 - 2/k$ Cut Partition | [ch04_multiway_kcut.cpp](approx_algorithms_cpp/ch04_multiway_kcut.cpp) | [ch04_multiway_kcut.py](approx_algorithms_python/ch04_multiway_kcut.py) |
| **5** | $k$-Center | 2-approx (Parametric Pruning) | [ch05_kcenter.cpp](approx_algorithms_cpp/ch05_kcenter.cpp) | [ch05_kcenter.py](approx_algorithms_python/ch05_kcenter.py) |
| **6** | Feedback Vertex Set | 2-approx (Local Ratio) | [ch06_feedback_vertex_set.cpp](approx_algorithms_cpp/ch06_feedback_vertex_set.cpp) | [ch06_feedback_vertex_set.py](approx_algorithms_python/ch06_feedback_vertex_set.py) |
| **7** | Shortest Superstring | 3-approx Cycle Cover | [ch07_shortest_superstring.cpp](approx_algorithms_cpp/ch07_shortest_superstring.cpp) | [ch07_shortest_superstring.py](approx_algorithms_python/ch07_shortest_superstring.py) |
| **8** | Knapsack | FPTAS (Value Rounding) | [ch08_knapsack.cpp](approx_algorithms_cpp/ch08_knapsack.cpp) | [ch08_knapsack.py](approx_algorithms_python/ch08_knapsack.py) |
| **9** | Bin Packing | Asymptotic PTAS (APTAS) | [ch09_bin_packing.cpp](approx_algorithms_cpp/ch09_bin_packing.cpp) | [ch09_bin_packing.py](approx_algorithms_python/ch09_bin_packing.py) |
| **10** | Minimum Makespan | PTAS (List Scheduling) | [ch10_makespan.cpp](approx_algorithms_cpp/ch10_makespan.cpp) | [ch10_makespan.py](approx_algorithms_python/ch10_makespan.py) |
| **11** | Euclidean TSP | PTAS (Quadtree dissection) | [ch11_euclidean_tsp.cpp](approx_algorithms_cpp/ch11_euclidean_tsp.cpp) | [ch11_euclidean_tsp.py](approx_algorithms_python/ch11_euclidean_tsp.py) |
| **12-14** | LP algorithms | $f$-approx Rounding & Dual Fitting | [ch12_14_lp_algorithms.cpp](approx_algorithms_cpp/ch12_14_lp_algorithms.cpp) | [ch12_14_lp_algorithms.py](approx_algorithms_python/ch12_14_lp_algorithms.py) |
| **13** | Knapsack FPTAS | FPTAS $(1-\epsilon)\cdot\text{OPT}$ | [ch13_knapsack.cpp](approx_algorithms_cpp/ch13_knapsack.cpp) | [ch13_knapsack.py](approx_algorithms_python/ch13_knapsack.py) |
| **15** | Weighted Vertex Cover | 2-approx (Primal-Dual Schema) | [ch15_weighted_vertex_cover_pd.cpp](approx_algorithms_cpp/ch15_weighted_vertex_cover_pd.cpp) | [ch15_weighted_vertex_cover_pd.py](approx_algorithms_python/ch15_weighted_vertex_cover_pd.py) |
| **16** | Randomized Rounding | $(1 - 1/e)$-approx Max-SAT | [ch16_randomized_rounding.cpp](approx_algorithms_cpp/ch16_randomized_rounding.cpp) | [ch16_randomized_rounding.py](approx_algorithms_python/ch16_randomized_rounding.py) |
| **17** | Chernoff Bounds | $O(\ln n)$ Set Cover rounding | [ch17_chernoff_bounds.cpp](approx_algorithms_cpp/ch17_chernoff_bounds.cpp) | [ch17_chernoff_bounds.py](approx_algorithms_python/ch17_chernoff_bounds.py) |
| **18** | SDP for Max-Cut | $\ge 0.878$-approx (Goemans-Williamson) | [ch18_sdp_maxcut.cpp](approx_algorithms_cpp/ch18_sdp_maxcut.cpp) | [ch18_sdp_maxcut.py](approx_algorithms_python/ch18_sdp_maxcut.py) |
| **19** | Multiway Cut LP Rounding | 1.3438-approx (CKR Rounding) | [ch19_multiway_rounding.cpp](approx_algorithms_cpp/ch19_multiway_rounding.cpp) | [ch19_multiway_rounding.py](approx_algorithms_python/ch19_multiway_rounding.py) |
| **21** | Steiner Forest | 2-approx (AKR Primal-Dual) | [ch21_steiner_forest.cpp](approx_algorithms_cpp/ch21_steiner_forest.cpp) | [ch21_steiner_forest.py](approx_algorithms_python/ch21_steiner_forest.py) |
| **22** | Steiner Network | 2-approx (Jain's Iterative Rounding) | [ch22_steiner_network.cpp](approx_algorithms_cpp/ch22_steiner_network.cpp) | [ch22_steiner_network.py](approx_algorithms_python/ch22_steiner_network.py) |
| **23** | Feedback Vertex Set | 2-approx (Primal-Dual Schema) | [ch23_primal_dual_fvs.cpp](approx_algorithms_cpp/ch23_primal_dual_fvs.cpp) | [ch23_primal_dual_fvs.py](approx_algorithms_python/ch23_primal_dual_fvs.py) |
| **24** | Facility Location | 3-approx Primal-Dual & Greedy | [ch24_facility_location.cpp](approx_algorithms_cpp/ch24_facility_location.cpp) | [ch24_facility_location.py](approx_algorithms_python/ch24_facility_location.py) |
| **26** | SDP for Max 2-SAT | $\ge 0.878$-approx (Goemans-Williamson) | [ch26_sdp_sat.cpp](approx_algorithms_cpp/ch26_sdp_sat.cpp) | [ch26_sdp_sat.py](approx_algorithms_python/ch26_sdp_sat.py) |
| **30** | Multicut in Trees | 2-approx LCA Depth | [ch30_tree_multicut.cpp](approx_algorithms_cpp/ch30_tree_multicut.cpp) | [ch30_tree_multicut.py](approx_algorithms_python/ch30_tree_multicut.py) |

*Note: Chapters 20, 25, 27, 28, and 29 focus entirely on complexity bounds and hardness proofs (such as the PCP theorem and Unique Games Conjecture) and thus do not contain programming implementations.*

---

## 💻 C++ Code Usage Example

To call the static library algorithms inside your own project, include the central header interface `approx_algorithms.hpp`:

```cpp
#include "approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;

    // Define a path graph P_4: 0-1-2-3
    Graph graph(4);
    graph[0] = {1};
    graph[1] = {0, 2};
    graph[2] = {1, 3};
    graph[3] = {2};

    // Run the factor-2 maximal matching-based vertex cover
    std::set<int> cover = vertex_cover_approx_2(graph);

    std::cout << "Cover vertices: ";
    for (int v : cover) {
        std::cout << v << " ";
    }
    std::cout << "\n"; // Outputs: 1 2
    return 0;
}
```

For more detailed code patterns (such as solving SDP or LP relaxations), inspect the files under the **[approx_algorithms_cpp/examples/](approx_algorithms_cpp/examples/)** directory.

---

## 📜 License & Citation

This project is licensed under the MIT License - see the LICENSE file for details.

If you use this companion codebase for teaching or research, please cite:
```text
Vazirani, Vijay V. "Approximation Algorithms". Springer-Verlag Berlin Heidelberg, 2001.
```

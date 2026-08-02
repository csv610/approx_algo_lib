# Approximation Algorithms: Companion Code for Vazirani's Book

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/compiler_support/23)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains complete, runnable implementation companions for Vijay Vazirani's seminal textbook, **"Approximation Algorithms" (Springer, 2001)**. 

It provides dual codebases (Python & C++23) to serve both educational and software engineering purposes, mapping mathematical proofs directly to executable computer code. The C++23 implementations are package-ready and backed by **Google OR-Tools** for industrial-strength Linear Programming.

---

## Key Features

* **Dual Implementations**:
  * **Python Codebase**: Highly readable, serving as interactive "executable pseudocode".
  * **C++23 Static Library (`aal`)**: A type-safe, modular implementation under `namespace aal` utilizing modern standard features (`std::print`, ranges, structured bindings).
* **Production-Grade Solvers**:
  - The C++ library integrates **Google OR-Tools (GLOP)** for solving LP relaxations efficiently (replacing standard textbook simplex implementations).
  - Includes a custom spherical gradient-projection Semidefinite Programming (SDP) optimizer built from scratch.
* **Zero Compilation Issues**: Highly portable, zero-warning builds with CMake.
* **100% Chapter Parity**: 30 companion chapters containing exact algorithms, runnable Python listings, and standalone C++ example binaries.

---

## Repository Directory Layout

```text
ApproxAlgo/
├── approx_algo.tex                  # LaTeX master for the companion book
├── chapters/*.tex                   # 30 chapter files, one per Vazirani chapter
├── backmatter/                      # Appendices (e.g. how-to-run)
├── src/
│   ├── pycodes/                     # Python Codebase (49 modules)
│   │   ├── main.py                  # Test suite runner: runs every chapter demo
│   │   ├── book.py                  # Chapter-to-module registry
│   │   └── *.py                     # One module per chapter + helper libraries
│   └── cppcodes/                     # C++23 Codebase (31 modules)
│       ├── CMakeLists.txt            # Build configuration
│       ├── approx_algorithms.hpp     # Unified public API header interface
│       ├── simplex.hpp               # Google OR-Tools GLOP solver wrapper
│       ├── main.cpp                  # Main test suite runner for C++
│       ├── USAGE.md                  # Integration guide & code reference
│       ├── examples/                 # 24 Standalone, compile-ready example files
│       │   ├── CMakeLists.txt
│       │   ├── ex_01_vertex_cover.cpp
│       │   └── ... ex_30_tree_multicut.cpp
│       └── tests/                    # Unit tests (GoogleTest)
├── README.md                         # Main repository overview
└── .gitignore                        # Residual compile & environment exclusions
```

---

## Quick Start & Build Guide

### Running Python Demos
Run all Python chapter demos synchronously:
```bash
cd src/pycodes
python3 main.py
```

### Building the C++ Library & Examples

#### Prerequisite: Google OR-Tools
Ensure Google OR-Tools is installed on your system:
* **macOS (Homebrew)**: `brew install or-tools`
* **Ubuntu/Debian**: `sudo apt-get install libortools-dev`

#### Building with CMake
```bash
cd src/cppcodes
cmake -B build -S .
cmake --build build
```
This compiles and creates:
* **`libapprox_algo.a`** (Static Library target `approx_algo`)
* **`approx_algo_demos`** (Global test runner executable `./build/approx_algo_demos`)
* **24 Standalone Example Binaries** (under `./build/examples/`)

### Building the Companion Book (LaTeX)
The companion notes compile with any standard TeX engine (e.g. `latexmk`):
```bash
latexmk -pdf approx_algo.tex
```
The book contains **30 chapters**, each with a self-contained, runnable Python
listing (dependencies: Python standard library only), plus a how-to-run appendix
in `backmatter/running_the_code.tex`.

---

## Implemented Chapters & Approximation Bounds

The repository's book covers all 30 chapters of Vazirani's book; the algorithmic chapters with runnable implementations are:

| Chapter | Algorithm / Topic | Approximation Bound | C++ Source | Python Source |
| :---: | :--- | :---: | :---: | :---: |
| **1** | Vertex Cover | Factor-2 (Maximal Matching) | [intro.cpp](src/cppcodes/intro.cpp) | [intro.py](src/pycodes/intro.py) |
| **2** | Set Cover | Greedy $H_n$-approx | [set_cover.cpp](src/cppcodes/set_cover.cpp) | [set_cover.py](src/pycodes/set_cover.py) |
| **3** | Steiner Tree & TSP | 1.5-approx (Christofides) | [steiner_tsp.cpp](src/cppcodes/steiner_tsp.cpp) | [steiner_tsp.py](src/pycodes/steiner_tsp.py) |
| **4** | Multiway Cut & $k$-Cut | $2 - 2/k$ Cut Partition | [multiway_kcut.cpp](src/cppcodes/multiway_kcut.cpp) | [multiway_kcut.py](src/pycodes/multiway_kcut.py) |
| **5** | $k$-Center | 2-approx (Parametric Pruning) | [kcenter.cpp](src/cppcodes/kcenter.cpp) | [kcenter.py](src/pycodes/kcenter.py) |
| **6** | Feedback Vertex Set | 2-approx (Local Ratio) | [feedback_vertex_set.cpp](src/cppcodes/feedback_vertex_set.cpp) | [feedback_vertex_set.py](src/pycodes/feedback_vertex_set.py) |
| **7** | Shortest Superstring | 3-approx Cycle Cover | [shortest_superstring.cpp](src/cppcodes/shortest_superstring.cpp) | [shortest_superstring.py](src/pycodes/shortest_superstring.py) |
| **8** | Knapsack | FPTAS (Value Rounding) | [knapsack.cpp](src/cppcodes/knapsack.cpp) | [knapsack.py](src/pycodes/knapsack.py) |
| **9** | Bin Packing | Asymptotic PTAS (APTAS) | [bin_packing.cpp](src/cppcodes/bin_packing.cpp) | [bin_packing.py](src/pycodes/bin_packing.py) |
| **10** | Minimum Makespan | PTAS (List Scheduling) | [makespan.cpp](src/cppcodes/makespan.cpp) | [makespan.py](src/pycodes/makespan.py) |
| **11** | Euclidean TSP | PTAS (Quadtree dissection) | [euclidean_tsp.cpp](src/cppcodes/euclidean_tsp.cpp) | [euclidean_tsp.py](src/pycodes/euclidean_tsp.py) |
| **12-14** | LP algorithms | $f$-approx Rounding & Dual Fitting | [lp_algorithms.cpp](src/cppcodes/lp_algorithms.cpp) | [lp_algorithms.py](src/pycodes/lp_algorithms.py) |
| **13** | Knapsack FPTAS | FPTAS $(1-\epsilon)\cdot\text{OPT}$ | [knapsack_ch13.cpp](src/cppcodes/knapsack_ch13.cpp) | [knapsack_ch13.py](src/pycodes/knapsack_ch13.py) |
| **15** | Weighted Vertex Cover | 2-approx (Primal-Dual Schema) | [weighted_vertex_cover_pd.cpp](src/cppcodes/weighted_vertex_cover_pd.cpp) | [weighted_vertex_cover_pd.py](src/pycodes/weighted_vertex_cover_pd.py) |
| **16** | Randomized Rounding | $(1 - 1/e)$-approx Max-SAT | [randomized_rounding.cpp](src/cppcodes/randomized_rounding.cpp) | [randomized_rounding.py](src/pycodes/randomized_rounding.py) |
| **17** | Chernoff Bounds | $O(\ln n)$ Set Cover rounding | [chernoff_bounds.cpp](src/cppcodes/chernoff_bounds.cpp) | [chernoff_bounds.py](src/pycodes/chernoff_bounds.py) |
| **18** | SDP for Max-Cut | $\ge 0.878$-approx (Goemans-Williamson) | [sdp_maxcut.cpp](src/cppcodes/sdp_maxcut.cpp) | [sdp_maxcut.py](src/pycodes/sdp_maxcut.py) |
| **19** | Multiway Cut LP Rounding | 1.3438-approx (CKR Rounding) | [multiway_rounding.cpp](src/cppcodes/multiway_rounding.cpp) | [multiway_rounding.py](src/pycodes/multiway_rounding.py) |
| **21** | Steiner Forest | 2-approx (AKR Primal-Dual) | [steiner_forest.cpp](src/cppcodes/steiner_forest.cpp) | [steiner_forest.py](src/pycodes/steiner_forest.py) |
| **22** | Steiner Network | 2-approx (Jain's Iterative Rounding) | [steiner_network.cpp](src/cppcodes/steiner_network.cpp) | [steiner_network.py](src/pycodes/steiner_network.py) |
| **23** | Feedback Vertex Set | 2-approx (Primal-Dual Schema) | [primal_dual_fvs.cpp](src/cppcodes/primal_dual_fvs.cpp) | [primal_dual_fvs.py](src/pycodes/primal_dual_fvs.py) |
| **24** | Facility Location | 3-approx Primal-Dual & Greedy | [facility_location.cpp](src/cppcodes/facility_location.cpp) | [facility_location.py](src/pycodes/facility_location.py) |
| **26** | SDP for Max 2-SAT | $\ge 0.878$-approx (Goemans-Williamson) | [sdp_sat.cpp](src/cppcodes/sdp_sat.cpp) | [sdp_sat.py](src/pycodes/sdp_sat.py) |
| **30** | Multicut in Trees | 2-approx LCA Depth | [tree_multicut.cpp](src/cppcodes/tree_multicut.cpp) | [tree_multicut.py](src/pycodes/tree_multicut.py) |

*Note: The companion book covers all 30 chapters of Vazirani. The table above lists the algorithmic chapters with full dual (Python + C++) runnable implementations; the remaining chapters center on complexity and hardness proofs (PCP theorem, quantifier-free inapproximability, Unique Games Conjecture) and, wherever an algorithmic core exists, are accompanied by a runnable Python listing.*

---

## C++ Code Usage Example

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

For more detailed code patterns (such as solving SDP or LP relaxations), inspect the files under the **[src/cppcodes/examples/](src/cppcodes/examples/)** directory.

---

## License & Citation

This project is licensed under the MIT License - see the LICENSE file for details.

If you use this companion codebase for teaching or research, please cite:
```text
Vazirani, Vijay V. "Approximation Algorithms". Springer-Verlag Berlin Heidelberg, 2001.
```

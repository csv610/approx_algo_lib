# How to Use the `aal` Library

The **`aal` (Approximate Algorithm Library)** is a modern C++23 static library implementing the companion algorithms for Vazirani's *Approximation Algorithms* book, backed by **Google OR-Tools** for industrial-strength Linear Programming.

---

## Prerequisite: Google OR-Tools
Before building `aal`, you must have Google OR-Tools installed:
* **macOS (Homebrew)**:
  ```bash
  brew install or-tools
  ```
* **Ubuntu/Debian**:
  ```bash
  sudo apt-get install libortools-dev
  ```

---

## 1. How to Compile the Library

You can compile the static library and run the demo suite using the provided CMake configuration.

### Option A: Standard Build (CMake)
```bash
cd approx_algorithms_cpp
cmake -B build -S .
cmake --build build
```
This builds:
* **`libapprox_algo.a`** (the static library target `approx_algo`)
* **`approx_algo_demos`** (the executable test runner `./build/approx_algo_demos`)

### Option B: Direct CLI Compilation
You can compile directly linking against Google OR-Tools:
```bash
clang++ -std=c++23 *.cpp -O3 -lortools -o my_demo_runner
```

---

## 2. Integrating `aal` into Your Own CMake Project

To use `aal` in your own project, add the `approx_algorithms_cpp` directory as a subdirectory:

1. In your `CMakeLists.txt`:
   ```cmake
   # Add the library directory
   add_subdirectory(path/to/approx_algorithms_cpp)

   # Create your executable
   add_executable(my_program main.cpp)

   # Link against the approx_algo static library
   target_link_libraries(my_program PRIVATE approx_algo)
   ```

2. In your source code (`main.cpp`), include the public interface:
   ```cpp
   #include "approx_algorithms.hpp"
   ```

---

## 3. Basic Code Examples

Here are some code snippets showing how to call the primary algorithms. All types and functions are defined under the namespace **`aal`**.

### Example 1: Vertex Cover (Factor-2 Approximation)
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

### Example 2: Knapsack FPTAS
```cpp
#include "approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;

    std::vector<int> weights = {10, 20, 30, 40, 50};
    std::vector<int> values = {60, 100, 120, 200, 250};
    int capacity = 100;
    double epsilon = 0.1; // 10% tolerance from OPT

    // Run Knapsack FPTAS
    auto [selected_items, value] = knapsack_fptas(weights, values, capacity, epsilon);

    std::cout << "Selected Knapsack Value: " << value << "\n";
    std::cout << "Selected Items: ";
    for (int idx : selected_items) {
        std::cout << idx << " ";
    }
    std::cout << "\n";
    return 0;
}
```

### Example 3: Metric TSP (3/2 Christofides Heuristic)
```cpp
#include "approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;

    // Complete Graph K_4 with metric edge weights
    int n = 4;
    std::vector<Edge> edges = {{0,1}, {1,2}, {2,3}, {3,0}, {0,2}, {1,3}};
    std::vector<double> costs = {2.0, 2.0, 2.0, 2.0, 3.0, 3.0};

    // Run Christofides 3/2-approximation for Metric TSP
    std::vector<int> tour = tsp_christofides_1_5_approx(n, edges, costs);

    std::cout << "TSP Tour: ";
    for (int city : tour) {
        std::cout << city << " -> ";
    }
    std::cout << tour[0] << "\n";
    return 0;
}
```

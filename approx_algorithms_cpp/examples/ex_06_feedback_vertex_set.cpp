/**
 * Chapter 6: Feedback Vertex Set (2-Approximation)
 * 
 * Theory:
 *   Given an undirected graph with vertex weights, we seek a minimum weight set
 *   of vertices whose removal makes the graph a forest (acyclic).
 *   This implementation uses the local ratio / primal-dual approach to obtain
 *   a factor-2 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 6: Feedback Vertex Set (2-approx) ===\n";

    Multigraph g = {{0, {1, 2}}, {1, {0, 2}}, {2, {0, 1}}};
    std::map<int, double> weights = {{0, 2.0}, {1, 3.0}, {2, 1.5}};

    std::cout << "  Input Graph: Triangle 0-1-2 (contains one cycle)\n";

    auto fvs = feedback_vertex_set_approx(g, weights);

    std::cout << "  Selected Feedback Vertices to Remove: [";
    bool first = true;
    for (int v : fvs) {
        std::cout << (first ? "" : ", ") << v;
        first = false;
    }
    std::cout << "]\n\n";
    return 0;
}

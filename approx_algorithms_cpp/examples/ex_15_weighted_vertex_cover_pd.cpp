/**
 * Chapter 15: Weighted Vertex Cover via Primal-Dual
 * 
 * Theory:
 *   For graphs with vertex weights w_v, the primal-dual algorithm grows dual variables
 *   y_e associated with edges until a vertex constraint becomes tight. The tight vertices
 *   form a vertex cover, guaranteeing a factor-2 approximation.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 15: Weighted Vertex Cover (Primal-Dual) ===\n";

    std::vector<int> vertices = {0, 1, 2};
    std::vector<Edge> edges = {{0, 1}, {1, 2}};
    std::map<int, double> weights = {{0, 1.5}, {1, 2.0}, {2, 1.0}};

    auto [cover, duals] = vertex_cover_primal_dual(vertices, edges, weights);

    std::cout << "  Vertices Cover Size: " << cover.size() << "\n\n";
    return 0;
}

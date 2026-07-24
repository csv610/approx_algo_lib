"""
Chapter 18: Semidefinite Programming (Max-Cut)
===============================================
Vazirani Ch. 18: Goemans-Williamson Randomized Rounding for Max-Cut.
Implements:
1. Spherical gradient projection optimizer to find unit vector embeddings.
2. Randomized hyperplane rounding to partition vertices.
3. Verification of the 0.878 approximation ratio.
"""

import math
import random
from typing import List, Tuple, Dict

def compute_cut_weight(
    assignment: List[int],  # 0 or 1 for each vertex
    edges: List[Tuple[int, int]],
    weights: List[float]
) -> float:
    """Calculate the total weight of edges crossing the cut."""
    weight = 0.0
    for (u, v), w in zip(edges, weights):
        if assignment[u] != assignment[v]:
            weight += w
    return weight

def optimize_max_cut_vectors(
    n: int,
    edges: List[Tuple[int, int]],
    weights: List[float],
    dim: int = 8,
    lr: float = 0.05,
    epochs: int = 200
) -> List[List[float]]:
    """
    Finds unit vector embeddings in R^dim for Max-Cut.
    Uses gradient descent with projection to maximize:
      F = 1/2 * sum_{(u,v) in E} w_{uv} * (1 - v_u . v_v)
      
    This is equivalent to minimizing sum_{(u,v) in E} w_{uv} * (v_u . v_v)
    subject to ||v_i|| = 1.
    """
    # Initialize random unit vectors in R^dim
    v = []
    for _ in range(n):
        vec = [random.gauss(0.0, 1.0) for _ in range(dim)]
        mag = math.sqrt(sum(x*x for x in vec))
        v.append([x/mag for x in vec])
        
    # Gradient descent with spherical projection
    for _ in range(epochs):
        new_v = []
        for i in range(n):
            grad = [0.0] * dim
            # Compute gradient dF/dv_i = sum_{j: (i,j) in E} w_{ij} * v_j
            for (u, val_v), weight in zip(edges, weights):
                if u == i:
                    for k in range(dim):
                        grad[k] += weight * v[val_v][k]
                elif val_v == i:
                    for k in range(dim):
                        grad[k] += weight * v[u][k]
                        
            # Update: we want to minimize dot products, so we subtract gradient
            updated = [v[i][k] - lr * grad[k] for k in range(dim)]
            mag = math.sqrt(sum(x*x for x in updated))
            if mag < 1e-9:
                new_v.append(v[i]) # keep old if collapsed
            else:
                new_v.append([x/mag for x in updated])
        v = new_v
    return v

def goemans_williamson_max_cut(
    n: int,
    edges: List[Tuple[int, int]],
    weights: List[float],
    vectors: List[List[float]],
    trials: int = 500
) -> Tuple[List[int], float, float]:
    """
    Applies randomized hyperplane rounding (Goemans-Williamson).
    Chooses a random unit vector r in R^dim and sets assignment[i] based on sign of v_i . r.
    
    Returns:
        best_assignment: List of 0/1 values for each vertex.
        avg_weight: Average weight of cuts generated.
        best_weight: Maximum weight of cut generated.
    """
    dim = len(vectors[0])
    best_weight = -1.0
    best_assignment = []
    total_weight = 0.0
    
    for _ in range(trials):
        # Generate random unit vector r in R^dim
        r = [random.gauss(0.0, 1.0) for _ in range(dim)]
        mag = math.sqrt(sum(x*x for x in r))
        r = [x/mag for x in r]
        
        # Rounded assignment: 1 if dot product >= 0, else 0
        assignment = []
        for i in range(n):
            dot = sum(vectors[i][k] * r[k] for k in range(dim))
            assignment.append(1 if dot >= 0.0 else 0)
            
        weight = compute_cut_weight(assignment, edges, weights)
        total_weight += weight
        if weight > best_weight:
            best_weight = weight
            best_assignment = assignment
            
    return best_assignment, total_weight / trials, best_weight

def demo_sdp_max_cut():
    print("=" * 60)
    print("Chapter 18: Semidefinite Programming (Max-Cut)")
    print("=" * 60)
    
    # Example 1: 5-cycle C_5 (classic non-bipartite graph)
    # Optimal Max-Cut is 4
    n1 = 5
    edges1 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
    weights1 = [1.0, 1.0, 1.0, 1.0, 1.0]
    
    print("\n1. 5-Cycle Graph C_5 (unit weights):")
    print(f"  Vertices: {list(range(n1))}")
    print(f"  Edges:    {edges1}")
    
    vectors1 = optimize_max_cut_vectors(n1, edges1, weights1, dim=2, lr=0.1, epochs=300)
    
    # Calculate angles between adjacent nodes on the 2D plane
    print("\n  Optimized 2D Vector Embeddings on Circle:")
    for i in range(n1):
        print(f"    v_{i}: {[round(x, 4) for x in vectors1[i]]}")
        
    best_cut, avg_cut, max_cut = goemans_williamson_max_cut(n1, edges1, weights1, vectors1, trials=1000)
    
    print(f"\n  Exact Optimal Max-Cut Value:  4.0")
    print(f"  GW Rounding Average Cut Value: {avg_cut:.4f} (ratio to OPT: {avg_cut/4.0:.4f}, bound: 0.878)")
    print(f"  GW Rounding Best Cut Value:    {max_cut:.4f} (assignment: {best_cut})")
    
    # Example 2: Petersen Graph (10 vertices, 15 edges)
    # Optimal Max-Cut is 12
    n2 = 10
    edges2 = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),  # outer cycle
        (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),  # inner star
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)   # spokes
    ]
    weights2 = [1.0] * len(edges2)
    
    print("\n2. Petersen Graph (10 vertices, 15 edges):")
    vectors2 = optimize_max_cut_vectors(n2, edges2, weights2, dim=8, lr=0.05, epochs=300)
    best_cut2, avg_cut2, max_cut2 = goemans_williamson_max_cut(n2, edges2, weights2, vectors2, trials=1000)
    
    print(f"  Exact Optimal Max-Cut Value:  12.0")
    print(f"  GW Rounding Average Cut Value: {avg_cut2:.4f} (ratio to OPT: {avg_cut2/12.0:.4f}, bound: 0.878)")
    print(f"  GW Rounding Best Cut Value:    {max_cut2:.4f} (assignment: {best_cut2})")

if __name__ == "__main__":
    demo_sdp_max_cut()

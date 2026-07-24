"""
Chapter 15: Weighted Vertex Cover via Primal-Dual
===================================================
Vazirani Ch. 15: Weighted Vertex Cover using the Primal-Dual Schema.
Implements:
1. Primal-Dual 2-approximation algorithm.
2. Verification of primal cost <= 2 * dual objective.
"""

from typing import List, Tuple, Dict, Set

def vertex_cover_primal_dual(
    vertices: List[int],
    edges: List[Tuple[int, int]],
    weights: Dict[int, float]
) -> Tuple[Set[int], Dict[Tuple[int, int], float]]:
    """
    Primal-Dual 2-approximation algorithm for Weighted Vertex Cover.
    
    Returns:
        cover: Set of selected vertices.
        y: Dual variables assigned to each edge.
    """
    # Initialize dual variables y_e = 0 for each edge
    # To keep edge keys consistent, store them sorted (u, v) where u < v
    y = {}
    for u, v in edges:
        e = (min(u, v), max(u, v))
        y[e] = 0.0
        
    # Track the dual sum accumulated at each vertex: sum_{e: e \ni u} y_e
    vertex_dual_sums = {v: 0.0 for v in vertices}
    
    cover = set()
    
    # Process edges
    for u, v in edges:
        e = (min(u, v), max(u, v))
        
        # If neither endpoint is in the cover, we must cover this edge
        if u not in cover and v not in cover:
            # How much can we increase y_e?
            # It is constrained by the remaining weight of the two endpoints
            slack_u = weights[u] - vertex_dual_sums[u]
            slack_v = weights[v] - vertex_dual_sums[v]
            
            # Find the minimum slack
            raise_amount = min(slack_u, slack_v)
            
            # Increase dual variable
            y[e] += raise_amount
            vertex_dual_sums[u] += raise_amount
            vertex_dual_sums[v] += raise_amount
            
            # Add endpoints that become tight to the cover
            if abs(vertex_dual_sums[u] - weights[u]) < 1e-9:
                cover.add(u)
            if abs(vertex_dual_sums[v] - weights[v]) < 1e-9:
                cover.add(v)
                
    return cover, y

def demo_weighted_vertex_cover_pd():
    print("=" * 60)
    print("Chapter 15: Weighted Vertex Cover via Primal-Dual")
    print("=" * 60)
    
    # Example 1: Bipartite-like graph with customized weights
    vertices1 = [0, 1, 2, 3]
    edges1 = [(0, 1), (1, 2), (2, 3), (3, 0)]
    weights1 = {0: 3.0, 1: 2.0, 2: 4.0, 3: 1.5}
    
    print("\n1. Input Graph: 4-cycle C_4")
    print(f"  Weights: {weights1}")
    
    cover, y = vertex_cover_primal_dual(vertices1, edges1, weights1)
    
    # Calculate costs
    primal_cost = sum(weights1[v] for v in cover)
    dual_cost = sum(y.values())
    
    print(f"  Primal Cover:    {cover} (cost: {primal_cost:.2f})")
    print(f"  Dual Variables:  {y} (dual obj: {dual_cost:.2f})")
    print(f"  Primal/Dual Ratio: {primal_cost/dual_cost:.4f} (theoretical bound <= 2.00)")
    
    # Example 2: Star graph (center has high weight, leaves have low weight)
    # Center = 0, leaves = 1, 2, 3
    vertices2 = [0, 1, 2, 3]
    edges2 = [(0, 1), (0, 2), (0, 3)]
    weights2 = {0: 10.0, 1: 4.0, 2: 4.0, 3: 4.0}
    
    print("\n2. Star Graph (center weight 10, leaves weight 4)")
    cover2, y2 = vertex_cover_primal_dual(vertices2, edges2, weights2)
    primal_cost2 = sum(weights2[v] for v in cover2)
    dual_cost2 = sum(y2.values())
    
    print(f"  Primal Cover:    {cover2} (cost: {primal_cost2:.2f})")
    print(f"  Dual Variables:  {y2} (dual obj: {dual_cost2:.2f})")
    print(f"  Primal/Dual Ratio: {primal_cost2/dual_cost2:.4f}")

if __name__ == "__main__":
    demo_weighted_vertex_cover_pd()

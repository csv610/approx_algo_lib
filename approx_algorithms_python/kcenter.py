"""
Chapter 5: k-Center Problem
============================
Vazirani Ch. 5:
- Parametric pruning for metric k-center (2-approx)
- Weighted k-center (3-approx)
"""

from typing import Dict, List, Set, Tuple, Optional
import math


# Type aliases
Graph = Dict[int, Dict[int, float]]


def build_graphs_by_threshold(graph: Graph) -> Tuple[List[float], List[Graph]]:
    """
    Build graphs G_1, G_2, ..., G_m where G_i has edges with weight <= threshold_i.
    Returns sorted unique edge weights and corresponding graphs.
    """
    edges = []
    for u in graph:
        for v, w in graph[u].items():
            if u < v:
                edges.append(w)
    thresholds = sorted(set(edges))
    
    graphs = []
    for t in thresholds:
        g = {u: {} for u in graph}
        for u in graph:
            for v, w in graph[u].items():
                if w <= t:
                    g[u][v] = w
        graphs.append(g)
    
    return thresholds, graphs


def graph_square(graph: Graph) -> Graph:
    """Return G^2: edges between vertices at distance <= 2."""
    vertices = list(graph.keys())
    n = len(vertices)
    g2 = {u: {} for u in vertices}
    
    for u in vertices:
        # Direct neighbors
        for v in graph[u]:
            g2[u][v] = graph[u][v]
        # Distance 2 neighbors
        for v in graph[u]:
            for w in graph[v]:
                if w != u:
                    g2[u][w] = min(g2[u].get(w, float('inf')), 
                                    graph[u][v] + graph[v][w])
    
    return g2


def maximal_independent_set(graph: Graph) -> Set[int]:
    """Find maximal independent set greedily."""
    independent = set()
    remaining = set(graph.keys())
    
    while remaining:
        v = next(iter(remaining))
        independent.add(v)
        # Remove v and its neighbors
        to_remove = {v}
        to_remove.update(graph[v].keys())
        remaining -= to_remove
    
    return independent


def dominating_set_from_mis(graph: Graph, mis: Set[int]) -> Set[int]:
    """From MIS in G^2, construct dominating set in G."""
    # For each v in MIS, take its neighbors in original G
    result = set()
    for v in mis:
        result.update(graph[v].keys())
    return result


def kcenter_parametric_pruning(graph: Graph, k: int) -> Tuple[Set[int], float]:
    """
    2-approximation for Metric k-Center (Algorithm 5.3 in Vazirani).
    
    Algorithm:
    1. Sort edges by weight, build G_1, ..., G_m
    2. For each G_i, compute G_i^2 (square graph)
    3. Find maximal independent set M_i in G_i^2
    4. Check if M_i is a dominating set in G_i
    5. Return first M_i that dominates G_i
    
    Approximation factor: 2
    """
    if k >= len(graph):
        return set(graph.keys()), 0.0
    
    thresholds, graphs = build_graphs_by_threshold(graph)
    
    for t, g in zip(thresholds, graphs):
        g2 = graph_square(g)
        mis = maximal_independent_set(g2)
        
        if len(mis) <= k:
            # Found solution with radius t
            centers = set(mis)
            # Ensure exactly k centers
            if len(centers) < k:
                remaining = set(graph.keys()) - centers
                centers.update(list(remaining)[:k - len(centers)])
            centers = set(list(centers)[:k])
            return centers, t
    
    # Fallback: all vertices
    return set(graph.keys()), max(thresholds) if thresholds else 0.0


def kcenter_2approx_dominating_set(graph: Graph, k: int) -> Tuple[Set[int], float]:
    """
    Alternative 2-approx for k-center: binary search on radius + dominating set.
    More practical implementation.
    """
    thresholds, graphs = build_graphs_by_threshold(graph)
    
    def has_dominating_set_of_size(g: Graph, k: int) -> Tuple[bool, Set[int]]:
        """Check if g has dominating set of size <= k using greedy."""
        # Greedy: repeatedly pick vertex covering most uncovered
        uncovered = set(g.keys())
        centers = set()
        
        while uncovered and len(centers) < k:
            best = max(uncovered, key=lambda v: len(set(g[v].keys()) & uncovered) + 1)
            centers.add(best)
            uncovered -= {best}
            uncovered -= set(g[best].keys())
        
        return len(uncovered) == 0, centers
    
    # Binary search on threshold
    left, right = 0, len(thresholds) - 1
    best_centers = set(graph.keys())
    best_radius = thresholds[-1] if thresholds else 0.0
    
    while left <= right:
        mid = (left + right) // 2
        t = thresholds[mid]
        g = graphs[mid]
        
        ok, centers = has_dominating_set_of_size(g, k)
        if ok:
            best_centers = centers
            best_radius = t
            right = mid - 1
        else:
            left = mid + 1
    
    return best_centers, best_radius


def weighted_kcenter_3approx(graph: Graph, weights: Dict[int, float], W: float) -> Tuple[Set[int], float]:
    """
    3-approximation for Weighted k-Center (Algorithm 5.10 in Vazirani).
    
    Problem: Find S ⊆ V with weight(S) ≤ W minimizing max_v min_{u∈S} cost(u,v)
    
    Algorithm:
    1. Build G_1^2, G_2^2, ... as before
    2. For each G_i^2, find maximal independent set M_i
    3. For each v ∈ M_i, assign weight si(v) = weight of min-weight vertex in N[v] in G_i
    3. Check if total weight of S_i = {si(v) | v ∈ M_i} ≤ W
    4. Return S_j for minimum j satisfying weight constraint
    
    Approximation factor: 3
    """
    if W >= sum(weights.values()):
        return set(graph.keys()), 0.0
    
    thresholds, graphs = build_graphs_by_threshold(graph)
    
    for t, g in zip(thresholds, graphs):
        g2 = graph_square(g)
        mis = maximal_independent_set(g2)
        
        # For each v in MIS, find min-weight vertex in its closed neighborhood in g
        selected = set()
        total_weight = 0.0
        
        for v in mis:
            # Neighborhood of v in g
            nbd = {v}
            nbd.update(g[v].keys())
            
            # Only consider vertices that exist in weights
            nbd = {u for u in nbd if u in weights}
            if not nbd:
                continue
                
            # Pick minimum weight vertex
            best = min(nbd, key=lambda x: weights[x])
            selected.add(best)
            total_weight += weights[best]
        
        if total_weight <= W:
            return selected, t
    
    return set(graph.keys()), max(thresholds) if thresholds else 0.0


def demo_kcenter():
    print("=" * 60)
    print("Chapter 5: k-Center Problem")
    print("=" * 60)
    
    # Example: cities on a line
    print("\n1. Unweighted k-Center (Parametric Pruning - 2-approx)")
    graph = {
        0: {1: 1, 2: 3},
        1: {0: 1, 2: 2, 3: 4},
        2: {0: 3, 1: 2, 3: 1, 4: 3},
        3: {1: 4, 2: 1, 4: 2},
        4: {2: 3, 3: 2}
    }
    for k in [1, 2, 3]:
        centers, radius = kcenter_parametric_pruning(graph, k)
        print(f"  k={k}: centers={centers}, radius={radius}")
    
    # Tight example: wheel graph
    print("\n2. Tight Example for 2-approx (Wheel Graph)")
    # Center 0 connected to all others with weight 1
    # Outer cycle with weight 2
    n = 6
    graph = {i: {} for i in range(n + 1)}
    for i in range(1, n + 1):
        graph[0][i] = 1.0
        graph[i][0] = 1.0
    for i in range(1, n + 1):
        j = (i % n) + 1
        graph[i][j] = 2.0
        graph[j][i] = 2.0
    
    for k in [1, 2]:
        centers, radius = kcenter_parametric_pruning(graph, k)
        print(f"  k={k}: centers={centers}, radius={radius}")
        if k == 1:
            print(f"  Optimal: center=0, radius=1. Approx ratio={radius/1.0}")
    
    # Weighted k-center
    print("\n3. Weighted k-Center (3-approx)")
    weights = {0: 1, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3}
    W = 4.0
    centers, radius = weighted_kcenter_3approx(graph, weights, W)
    print(f"  Vertex weights: {weights}")
    print(f"  Budget W={W}")
    print(f"  Selected centers: {centers}, weight={sum(weights[c] for c in centers)}, radius={radius}")


if __name__ == "__main__":
    demo_kcenter()
"""
Chapter 1: Introduction - Vertex Cover
======================================
Vazirani Ch. 1: Factor-2 approximation for Vertex Cover via maximal matching.

Key idea: Maximal matching gives a lower bound on OPT.
Algorithm: Find maximal matching, output all matched vertices.
Approximation factor: 2
"""

from typing import List, Set, Tuple
import random


# Type aliases
Edge = Tuple[int, int]
Graph = List[List[int]]  # Adjacency list


def maximal_matching(graph: Graph) -> List[Edge]:
    """
    Greedy maximal matching: pick edges greedily, remove endpoints.
    Returns list of matched edges.
    """
    n = len(graph)
    matched = [False] * n
    matching = []
    
    for u in range(n):
        if not matched[u]:
            for v in graph[u]:
                if not matched[v]:
                    matching.append((u, v))
                    matched[u] = matched[v] = True
                    break
    return matching


def vertex_cover_approx_2(graph: Graph) -> Set[int]:
    """
    Factor-2 approximation for Vertex Cover (Algorithm 1.2 in Vazirani).
    
    Algorithm:
    1. Find a maximal matching M
    2. Output all vertices incident to edges in M
    
    Approximation guarantee: |C| <= 2 * OPT
    Proof: |M| <= OPT (matching is lower bound), |C| = 2|M|
    """
    matching = maximal_matching(graph)
    cover = set()
    for u, v in matching:
        cover.add(u)
        cover.add(v)
    return cover


def vertex_cover_approx_2_edge_weighted(edges: List[Edge], weights: dict) -> Set[int]:
    """
    Weighted version using LP rounding (Chapter 14 in Vazirani).
    For now, simple greedy: pick highest weight edge endpoints.
    """
    # Greedy heuristic for weighted: pick min weight vertex cover
    # This is just a greedy heuristic, not a true 2-approx for weighted
    # True 2-approx uses LP rounding (Ch. 14)
    n = max(max(u, v) for u, v in edges) + 1
    covered = [False] * len(edges)
    cover = set()
    edges_sorted = sorted(edges, key=lambda e: weights[e[0]] + weights[e[1]])
    
    for u, v in edges_sorted:
        if not covered[edges.index((u, v))]:
            # Pick the cheaper endpoint
            if weights[u] < weights[v]:
                cover.add(u)
            else:
                cover.add(v)
            # Mark covered edges
            for i, (u2, v2) in enumerate(edges):
                if u == u2 or u == v2 or v == u2 or v == v2:
                    covered[i] = True
    return cover


# Tight example: Complete bipartite K_{n,n}
def tight_example_k_n_n(n: int) -> Graph:
    """Complete bipartite K_{n,n} - tight example for factor 2."""
    graph = [[] for _ in range(2 * n)]
    for u in range(n):
        for v in range(n, 2 * n):
            graph[u].append(v)
            graph[v].append(u)
    return graph


def vertex_cover_exact_bruteforce(graph: Graph) -> Set[int]:
    """Exact vertex cover via brute force (for small graphs only)."""
    n = len(graph)
    best = set(range(n))
    
    for mask in range(1 << n):
        cover = {i for i in range(n) if mask & (1 << i)}
        valid = True
        for u in range(n):
            for v in graph[u]:
                if u not in cover and v not in cover:
                    valid = False
                    break
            if not valid:
                break
        if valid and len(cover) < len(best):
            best = cover
    return best


def demo_vertex_cover():
    print("=" * 60)
    print("Chapter 1: Vertex Cover - Factor 2 Approximation")
    print("=" * 60)
    
    # Example 1: Complete bipartite K_{4,4} (tight example)
    print("\n1. Tight Example: K_{4,4}")
    g = tight_example_k_n_n(4)
    cover = vertex_cover_approx_2(g)
    exact = vertex_cover_exact_bruteforce(g)
    print(f"  Graph: K_{{4,4}} (8 vertices, 16 edges)")
    print(f"  Approx cover size: {len(cover)}")
    print(f"  Optimal cover size: {len(exact)}")
    print(f"  Ratio: {len(cover) / len(exact):.2f}")
    print(f"  (Optimal picks one side: 4 vertices)")
    
    # Example 2: Random graph
    print("\n2. Random Graph (10 vertices, p=0.3)")
    n = 10
    p = 0.3
    random.seed(42)
    g = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < p:
                g[i].append(j)
                g[j].append(i)
    
    cover = vertex_cover_approx_2(g)
    exact = vertex_cover_exact_bruteforce(g)
    print(f"  Graph: {n} vertices, {sum(len(adj) for adj in g)//2} edges")
    print(f"  Approx cover size: {len(cover)}")
    print(f"  Optimal cover size: {len(exact)}")
    print(f"  Ratio: {len(cover) / len(exact):.2f}")
    
    # Example 3: Path graph P_n
    print("\n3. Path Graph P_5")
    g = [[] for _ in range(5)]
    for i in range(4):
        g[i].append(i+1)
        g[i+1].append(i)
    cover = vertex_cover_approx_2(g)
    exact = vertex_cover_exact_bruteforce(g)
    print(f"  Approx cover: {sorted(cover)}")
    print(f"  Optimal cover: {sorted(exact)}")
    print(f"  Ratio: {len(cover) / len(exact):.2f}")


if __name__ == "__main__":
    demo_vertex_cover()
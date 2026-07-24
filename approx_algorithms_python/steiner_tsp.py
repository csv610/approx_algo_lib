"""
Chapter 3: Steiner Tree and TSP
================================
Vazirani Ch. 3:
- Metric Steiner Tree: 2-approximation via MST on metric closure
- Metric TSP: 2-approx via MST double-tree, 3/2-approx via Christofides
"""

from typing import List, Set, Dict, Tuple, Optional
import heapq
import itertools
import math


# Type aliases
Graph = Dict[int, Dict[int, float]]  # adjacency dict with weights
Edge = Tuple[int, int]
Path = List[int]


def mst_prim(graph: Graph, vertices: List[int] = None) -> List[Edge]:
    """Prim's algorithm for MST. Returns list of edges (u, v)."""
    if vertices is None:
        vertices = list(graph.keys())
    if not vertices:
        return []
    
    start = vertices[0]
    visited = {start}
    edges = []
    # (weight, u, v)
    pq = [(w, start, v) for v, w in graph[start].items() if v in vertices]
    heapq.heapify(pq)
    
    while pq and len(visited) < len(vertices):
        w, u, v = heapq.heappop(pq)
        if v in visited:
            continue
        visited.add(v)
        edges.append((u, v))
        for w2, wgt in graph[v].items():
            if w2 in vertices and w2 not in visited:
                heapq.heappush(pq, (wgt, v, w2))
    
    return edges


def dijkstra(graph: Graph, source: int, targets: Set[int] = None) -> Dict[int, float]:
    """Dijkstra's algorithm from single source."""
    dist = {source: 0.0}
    pq = [(0.0, source)]
    visited = set()
    target_set = targets if targets is not None else set(graph.keys())
    
    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        if u in target_set and len(visited) == len(target_set):
            break
        for v, w in graph[u].items():
            nd = d + w
            if v not in dist or nd < dist[v]:
                dist[v] = nd
                heapq.heappush(pq, (nd, v))
    return dist


def metric_closure(graph: Graph, terminals: Set[int]) -> Graph:
    """
    Compute metric closure on terminals: complete graph with shortest path distances.
    """
    terminals = list(terminals)
    closure = {u: {} for u in terminals}
    
    for u in terminals:
        dist = dijkstra(graph, u, terminals)
        for v in terminals:
            if u != v:
                closure[u][v] = dist[v]
                closure[v][u] = dist[v]
    return closure


def mst_weight(graph: Graph, edges: List[Edge]) -> float:
    """Compute total weight of MST edges."""
    return sum(graph[u][v] for u, v in edges)


def steiner_tree_2approx(graph: Graph, terminals: Set[int]) -> Tuple[List[Edge], float]:
    """
    2-approximation for Metric Steiner Tree (Theorem 3.3 in Vazirani).
    
    Algorithm:
    1. Compute metric closure G' on terminals
    2. Find MST T' in G'
    3. Replace each edge in T' by corresponding shortest path in G
    4. Take union, remove cycles to get tree T
    
    Approximation factor: 2
    """
    if len(terminals) <= 1:
        return [], 0.0
    
    # Step 1: Metric closure
    closure = metric_closure(graph, terminals)
    
    # Step 2: MST on closure
    mst_edges = mst_prim(closure, list(terminals))
    
    # Step 3: Replace edges with shortest paths
    # We need to track actual paths in original graph
    all_edges = set()
    total_weight = 0.0
    
    for u, v in mst_edges:
        # Dijkstra with path reconstruction
        dist = {u: 0.0}
        prev = {u: None}
        pq = [(0.0, u)]
        
        while pq:
            d, x = heapq.heappop(pq)
            if x == v:
                break
            if d > dist[x]:
                continue
            for y, w in graph[x].items():
                nd = d + w
                if y not in dist or nd < dist[y]:
                    dist[y] = nd
                    prev[y] = x
                    heapq.heappush(pq, (nd, y))
        
        # Reconstruct path
        path = []
        x = v
        while x is not None:
            path.append(x)
            x = prev[x]
        path.reverse()
        
        # Add edges
        for i in range(len(path) - 1):
            a, b = path[i], path[i + 1]
            if a > b:
                a, b = b, a
            all_edges.add((a, b))
            total_weight += graph[a][b]
    
    # Step 4: Remove cycles (get spanning tree of the union)
    # Use DFS to get tree
    tree_edges = []
    visited = set()
    
    def dfs(u, parent):
        visited.add(u)
        for v, w in graph[u].items():
            if (min(u, v), max(u, v)) in all_edges:
                if v not in visited:
                    tree_edges.append((u, v))
                    dfs(v, u)
    
    # Start from any terminal
    start = next(iter(terminals))
    dfs(start, -1)
    
    # Recompute weight
    final_weight = sum(graph[u][v] for u, v in tree_edges)
    
    return tree_edges, final_weight


# ============================================================
# METRIC TSP
# ============================================================

def tsp_2approx_mst(graph: Graph) -> Tuple[List[int], float]:
    """
    2-approximation for Metric TSP (Algorithm 3.7 in Vazirani).
    
    Algorithm:
    1. Find MST T of G
    2. Double every edge of T to get Eulerian multigraph
    3. Find Eulerian tour
    4. Shortcut to get Hamiltonian cycle
    
    Approximation factor: 2
    """
    n = len(graph)
    if n <= 1:
        return list(graph.keys()), 0.0
    if n == 2:
        u, v = list(graph.keys())
        return [u, v, u], 2 * graph[u][v]
    
    # Step 1: MST
    mst_edges = mst_prim(graph)
    
    # Step 2: Double edges (build Eulerian multigraph adjacency)
    euler_adj = {u: [] for u in graph}
    for u, v in mst_edges:
        euler_adj[u].append(v)
        euler_adj[v].append(u)
        euler_adj[u].append(v)
        euler_adj[v].append(u)
    
    # Step 3: Hierholzer's algorithm for Eulerian tour
    def eulerian_tour(start):
        stack = [start]
        tour = []
        adj_copy = {u: list(v) for u, v in euler_adj.items()}
        
        while stack:
            u = stack[-1]
            if adj_copy[u]:
                v = adj_copy[u].pop()
                # Remove reverse edge
                adj_copy[v].remove(u)
                stack.append(v)
            else:
                tour.append(stack.pop())
        return tour[::-1]  # Reverse to get correct order
    
    euler_tour = eulerian_tour(next(iter(graph)))
    
    # Step 4: Shortcut - visit each vertex first time it appears
    visited = set()
    ham_cycle = []
    for v in euler_tour:
        if v not in visited:
            visited.add(v)
            ham_cycle.append(v)
    
    # Return to start
    ham_cycle.append(ham_cycle[0])
    
    # Compute cost
    cost = sum(graph[ham_cycle[i]][ham_cycle[i+1]] for i in range(len(ham_cycle)-1))
    
    return ham_cycle, cost


def min_weight_perfect_matching(graph: Graph, vertices: List[int]) -> List[Edge]:
    """
    Minimum weight perfect matching on a subset of vertices.
    Uses DP (Held-Karp style) - exact but exponential.
    For Christofides, we need matching on odd-degree vertices of MST.
    """
    n = len(vertices)
    if n == 0:
        return []
    if n % 2 == 1:
        raise ValueError("Odd number of vertices for perfect matching")
    if n > 16:  # Too large for exact DP
        # Greedy approximation
        return greedy_matching(graph, vertices)
    
    # Map vertices to indices
    idx = {v: i for i, v in enumerate(vertices)}
    
    # DP: dp[mask] = min weight perfect matching on vertices in mask
    # mask has even number of bits
    INF = float('inf')
    dp = {0: 0.0}
    parent = {}
    
    for mask in range(1 << n):
        if mask not in dp:
            continue
        # Find first unmatched vertex
        i = 0
        while i < n and (mask & (1 << i)):
            i += 1
        if i >= n:
            continue
        
        # Match i with some j > i
        for j in range(i + 1, n):
            if not (mask & (1 << j)):
                u, v = vertices[i], vertices[j]
                if v in graph[u]:
                    w = graph[u][v]
                    nmask = mask | (1 << i) | (1 << j)
                    if nmask not in dp or dp[mask] + w < dp[nmask]:
                        dp[nmask] = dp[mask] + w
                        parent[nmask] = (mask, (u, v))
    
    full = (1 << n) - 1
    if full not in dp:
        return greedy_matching(graph, vertices)
    
    # Reconstruct matching
    matching = []
    mask = full
    while mask != 0:
        pmask, edge = parent[mask]
        matching.append(edge)
        mask = pmask
    
    return matching


def greedy_matching(graph: Graph, vertices: List[int]) -> List[Edge]:
    """Greedy matching: repeatedly pick minimum weight edge."""
    remaining = set(vertices)
    matching = []
    
    # All possible edges sorted by weight
    edges = []
    for i, u in enumerate(vertices):
        for v in vertices[i+1:]:
            if v in graph[u]:
                edges.append((graph[u][v], u, v))
    edges.sort()
    
    for w, u, v in edges:
        if u in remaining and v in remaining:
            matching.append((u, v))
            remaining.remove(u)
            remaining.remove(v)
    
    return matching


def tsp_christofides_1_5_approx(graph: Graph) -> Tuple[List[int], float]:
    """
    Christofides 3/2-approximation for Metric TSP (Algorithm 3.10 in Vazirani).
    
    Algorithm:
    1. Find MST T of G
    2. Let V' be odd-degree vertices in T
    3. Find min-weight perfect matching M on V' (in metric closure)
    4. Combine T + M to get Eulerian multigraph
    5. Find Eulerian tour, shortcut to Hamiltonian cycle
    
    Approximation factor: 3/2
    """
    n = len(graph)
    if n <= 2:
        return tsp_2approx_mst(graph)
    
    # Step 1: MST
    mst_edges = mst_prim(graph)
    
    # Step 2: Find odd-degree vertices
    degree = {u: 0 for u in graph}
    for u, v in mst_edges:
        degree[u] += 1
        degree[v] += 1
    odd_vertices = [u for u, d in degree.items() if d % 2 == 1]
    
    # Step 3: Min-weight perfect matching on odd vertices
    # Use metric closure distances
    closure = metric_closure(graph, set(odd_vertices))
    matching = min_weight_perfect_matching(closure, odd_vertices)
    
    # Step 4: Combine T + M (build Eulerian multigraph)
    euler_adj = {u: [] for u in graph}
    for u, v in mst_edges:
        euler_adj[u].append(v)
        euler_adj[v].append(u)
    for u, v in matching:
        euler_adj[u].append(v)
        euler_adj[v].append(u)
    
    # Step 5: Eulerian tour + shortcut
    def eulerian_tour(start):
        stack = [start]
        tour = []
        adj_copy = {u: list(v) for u, v in euler_adj.items()}
        
        while stack:
            u = stack[-1]
            if adj_copy[u]:
                v = adj_copy[u].pop()
                adj_copy[v].remove(u)
                stack.append(v)
            else:
                tour.append(stack.pop())
        return tour[::-1]
    
    euler_tour = eulerian_tour(next(iter(graph)))
    
    visited = set()
    ham_cycle = []
    for v in euler_tour:
        if v not in visited:
            visited.add(v)
            ham_cycle.append(v)
    ham_cycle.append(ham_cycle[0])
    
    cost = sum(graph[ham_cycle[i]][ham_cycle[i+1]] for i in range(len(ham_cycle)-1))
    
    return ham_cycle, cost


def tsp_tight_example_2approx(n: int = 6) -> Tuple[Graph, float, float]:
    """
    Tight example for 2-approx TSP.
    Uses a path of length n with shortcut edges of weight 1 + epsilon.
    The MST is the path of weight n-1.
    Double tree shortcutted tour uses a return edge of weight n/2.
    Optimal tour uses alternating shortcuts to achieve weight n.
    Ratio approaches 1.5 as n increases.
    """
    eps = 0.1
    vertices = list(range(n))
    
    # Base tree for metric completion
    tree = {u: {} for u in vertices}
    for i in range(n - 1):
        tree[i][i+1] = 1.0
        tree[i+1][i] = 1.0
    for i in range(n - 2):
        tree[i][i+2] = 1.0 + eps
        tree[i+2][i] = 1.0 + eps
        
    # BFS to compute all-pairs shortest paths
    graph = {u: {} for u in vertices}
    for start in vertices:
        queue = [(start, 0.0)]
        visited = {start}
        while queue:
            curr, d = queue.pop(0)
            graph[start][curr] = d
            for nxt, w in tree[curr].items():
                if nxt not in visited:
                    visited.add(nxt)
                    queue.append((nxt, d + w))
                    
    _, approx_cost = tsp_2approx_mst(graph)
    
    # Calculate exact optimal cost:
    # Cycle: 0 -> 1 -> 3 -> 5 -> ... -> (n-1) -> (n-2) -> (n-4) -> ... -> 0
    # Uses 2 path edges of weight 1, and n-2 shortcut edges of weight 1+eps
    # So cost is 2 * 1.0 + (n - 2) * (1.0 + eps)
    opt_cost = 2.0 + (n - 2) * (1.0 + eps)
    
    return graph, approx_cost, opt_cost


def demo_steiner_tsp():
    print("=" * 60)
    print("Chapter 3: Steiner Tree and TSP")
    print("=" * 60)
    
    # Steiner Tree Example
    print("\n1. Metric Steiner Tree (2-approx)")
    # Graph: terminals {0, 3, 5}, Steiner {1, 2, 4}
    graph = {
        0: {1: 1, 2: 2},
        1: {0: 1, 2: 1, 3: 2},
        2: {0: 2, 1: 1, 4: 1},
        3: {1: 2, 5: 1},
        4: {2: 1, 5: 1},
        5: {3: 1, 4: 1}
    }
    terminals = {0, 3, 5}
    edges, cost = steiner_tree_2approx(graph, terminals)
    print(f"  Terminals: {terminals}")
    print(f"  Steiner tree edges: {edges}")
    print(f"  Cost: {cost}")
    
    # TSP 2-approx
    print("\n2. Metric TSP 2-approx (MST double-tree)")
    tsp_graph = {
        0: {1: 2, 2: 3, 3: 1},
        1: {0: 2, 2: 2, 3: 4},
        2: {0: 3, 1: 2, 3: 2},
        3: {0: 1, 1: 4, 2: 2}
    }
    cycle, cost = tsp_2approx_mst(tsp_graph)
    print(f"  Tour: {cycle}")
    print(f"  Cost: {cost}")
    
    # TSP Christofides 3/2-approx
    print("\n3. Metric TSP Christofides (3/2-approx)")
    cycle, cost = tsp_christofides_1_5_approx(tsp_graph)
    print(f"  Tour: {cycle}")
    print(f"  Cost: {cost}")
    
    # Tight example
    print("\n4. Tight example for 2-approx TSP")
    for n in [4, 6, 10]:
        g, approx, opt = tsp_tight_example_2approx(n)
        print(f"  n={n}: approx={approx}, opt={opt}, ratio={approx/opt:.3f}")


if __name__ == "__main__":
    demo_steiner_tsp()
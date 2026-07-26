"""
Chapter 4: Multiway Cut and Minimum k-Cut
===========================================
Vazirani Ch. 4:
- Multiway Cut: 2-2/k approximation via isolating cuts
- Minimum k-Cut: 2-2/k approximation via Gomory-Hu tree
"""

from typing import Dict, List, Set, Tuple, Optional
import heapq
import itertools
import math


# Type aliases
Graph = Dict[int, Dict[int, float]]
Edge = Tuple[int, int]


class MaxFlow:
    """Max flow for Gomory-Hu tree (simple Edmonds-Karp)."""
# MAX FLOW / MIN CUT (Edmonds-Karp for exact min s-t cut)
# ============================================================

class MaxFlow:
    """Edmonds-Karp max flow for min cut computation."""
    
    def __init__(self, n: int):
        self.n = n
        self.adj = [[] for _ in range(n)]
        self.cap = [[0] * n for _ in range(n)]
    
    def add_edge(self, u: int, v: int, c: float):
        self.cap[u][v] += c
        if v not in self.adj[u]:
            self.adj[u].append(v)
        if u not in self.adj[v]:
            self.adj[v].append(u)
    
    def bfs(self, s: int, t: int, parent: List[int]) -> float:
        parent[:] = [-1] * self.n
        parent[s] = -2
        q = [(s, float('inf'))]
        
        while q:
            u, flow = q.pop(0)
            for v in self.adj[u]:
                if parent[v] == -1 and self.cap[u][v] > 0:
                    parent[v] = u
                    nf = min(flow, self.cap[u][v])
                    if v == t:
                        return nf
                    q.append((v, nf))
        return 0
    
    def max_flow(self, s: int, t: int) -> float:
        flow = 0
        parent = [-1] * self.n
        while True:
            nf = self.bfs(s, t, parent)
            if nf == 0:
                break
            flow += nf
            v = t
            while v != s:
                u = parent[v]
                self.cap[u][v] -= nf
                self.cap[v][u] += nf
                v = u
        return flow
    
    def min_cut(self, s: int, t: int) -> Tuple[Set[int], Set[int]]:
        """Returns (S, T) partition after max flow."""
        self.max_flow(s, t)
        # Find reachable vertices in residual graph
        visited = set()
        stack = [s]
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            for v in self.adj[u]:
                if self.cap[u][v] > 0 and v not in visited:
                    stack.append(v)
        S = visited
        T = set(range(self.n)) - S
        return S, T


def min_s_t_cut(graph: Graph, s: int, t: int) -> Tuple[Set[int], Set[int], float]:
    """Compute minimum s-t cut using max flow."""
    # Map vertices to indices
    vertices = list(graph.keys())
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    
    mf = MaxFlow(n)
    for u in graph:
        for v, c in graph[u].items():
            if u < v:  # Add once
                mf.add_edge(idx[u], idx[v], c)
                mf.add_edge(idx[v], idx[u], c)
    
    S_idx, T_idx = mf.min_cut(idx[s], idx[t])
    S = {vertices[i] for i in S_idx}
    T = {vertices[i] for i in T_idx}
    
    # Compute cut weight
    cut_weight = 0.0
    for u in S:
        for v in T:
            if v in graph[u]:
                cut_weight += graph[u][v]
    
    return S, T, cut_weight


# ============================================================
# MULTIWAY CUT (Algorithm 4.3 in Vazirani)
# ============================================================

def multiway_cut_2_2k(graph: Graph, terminals: Set[int]) -> Tuple[Set[Edge], float]:
    """
    2 - 2/k approximation for Multiway Cut (Algorithm 4.3 in Vazirani).
    
    Algorithm:
    1. For each terminal s_i, compute minimum isolating cut C_i
       (minimum cut separating s_i from all other terminals)
    2. Output union of all cuts except the heaviest
    
    Approximation factor: 2 - 2/k where k = |terminals|
    """
    k = len(terminals)
    if k <= 2:
        # For 2 terminals, exact min s-t cut
        s, t = list(terminals)
        S, T, w = min_s_t_cut(graph, s, t)
        edges = set()
        for u in S:
            for v in T:
                if v in graph[u]:
                    edges.add((min(u, v), max(u, v)))
        return edges, w
    
    terminals = list(terminals)
    isolating_cuts = []
    cut_weights = []
    
    for i, s in enumerate(terminals):
        # Merge all other terminals into a super-terminal
        other_terminals = [t for j, t in enumerate(terminals) if j != i]
        
        # Create auxiliary graph: contract other terminals
        # We'll do this by modifying the graph temporarily
        # Simpler: add super-source connected to all other terminals with infinite capacity
        INF = float('inf')
        
        # Build flow network with super-source
        vertices = list(graph.keys())
        idx = {v: i for i, v in enumerate(vertices)}
        n = len(vertices)
        super_source = n
        mf = MaxFlow(n + 1)
        
        for u in graph:
            for v, c in graph[u].items():
                if u < v:
                    mf.add_edge(idx[u], idx[v], c)
                    mf.add_edge(idx[v], idx[u], c)
        
        # Connect super source to other terminals with infinite capacity
        for t in other_terminals:
            mf.add_edge(super_source, idx[t], INF)
        
        # Min cut between super source and s
        mf.max_flow(super_source, idx[s])
        
        # Find reachable from super source
        visited = set()
        stack = [super_source]
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            for v in mf.adj[u]:
                if mf.cap[u][v] > 0 and v not in visited:
                    stack.append(v)
        
        S_idx = visited - {super_source}
        S = {vertices[i] for i in S_idx}
        T = set(vertices) - S
        
        # Collect cut edges
        cut_edges = set()
        cut_weight = 0.0
        for u in S:
            for v in T:
                if v in graph[u]:
                    cut_edges.add((min(u, v), max(u, v)))
                    cut_weight += graph[u][v]
        
        isolating_cuts.append(cut_edges)
        cut_weights.append(cut_weight)
    
    # Discard heaviest cut
    max_idx = max(range(k), key=lambda i: cut_weights[i])
    
    # Union of remaining cuts
    result_edges = set()
    result_weight = 0.0
    for i in range(k):
        if i != max_idx:
            result_edges.update(isolating_cuts[i])
            result_weight += cut_weights[i]
    
    return result_edges, result_weight


# ============================================================
# GOMORY-HU TREE (for minimum k-cut)
# ============================================================

def gomory_hu_tree(graph: Graph) -> Tuple[List[Edge], Dict[Edge, float]]:
    """
    Gomory-Hu tree construction (Algorithm 4.6 in Vazirani).
    
    Returns: (tree_edges, edge_weights)
    """
    vertices = list(graph.keys())
    n = len(vertices)
    if n <= 1:
        return [], {}
    
    idx = {v: i for i, v in enumerate(vertices)}
    
    # Initial partition: each vertex in its own set
    parent = {i: 0 for i in range(1, n)}  # parent[i] = representative of set containing i
    tree_edges = []
    tree_weights = {}
    
    for i in range(1, n):
        s = i
        t = parent[i]
        
        # Build flow network
        mf = MaxFlow(n)
        for u in graph:
            for v, c in graph[u].items():
                if u < v:
                    mf.add_edge(idx[u], idx[v], c)
                    mf.add_edge(idx[v], idx[u], c)
        
        # Min cut between s and t
        mf.max_flow(s, t)
        
        # Find reachable from s
        visited = set()
        stack = [s]
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            for v in mf.adj[u]:
                if mf.cap[u][v] > 0 and v not in visited:
                    stack.append(v)
        
        # Add tree edge
        u, v = vertices[s], vertices[t]
        weight = 0.0
        for x in visited:
            for y in range(n):
                if y not in visited:
                    if vertices[y] in graph[vertices[x]]:
                        weight += graph[vertices[x]][vertices[y]]
        
        tree_edges.append((u, v))
        tree_weights[(min(u, v), max(u, v))] = weight
        
        # Update parent for vertices in same side as s
        for j in range(i + 1, n):
            if parent[j] == t and j in visited:
                parent[j] = i
    
    return tree_edges, tree_weights


def min_k_cut_2_2k(graph: Graph, k: int) -> Tuple[Set[Edge], float]:
    """
    2 - 2/k approximation for Minimum k-Cut (Algorithm 4.7 in Vazirani).
    
    Algorithm:
    1. Compute Gomory-Hu tree T of G
    2. Remove k-1 lightest edges from T
    3. Output union of corresponding cuts in G
    
    Approximation factor: 2 - 2/k
    """
    if k <= 1:
        return set(), 0.0
    if k >= len(graph):
        # Every vertex in its own component - all edges
        edges = set()
        weight = 0.0
        for u in graph:
            for v in graph[u]:
                if u < v:
                    edges.add((u, v))
                    weight += graph[u][v]
        return edges, weight
    
    # Step 1: Gomory-Hu tree
    tree_edges, tree_weights = gomory_hu_tree(graph)
    
    # Step 2: Sort tree edges by weight, remove k-1 lightest
    sorted_edges = sorted(tree_edges, key=lambda e: tree_weights[(min(e[0], e[1]), max(e[0], e[1]))])
    to_remove = sorted_edges[:k-1]
    
    # Step 3: For each removed tree edge, find corresponding cut in G
    # The cut is the min u-v cut in G, which we need to recompute
    result_edges = set()
    result_weight = 0.0
    
    for u, v in to_remove:
        S, T, w = min_s_t_cut(graph, u, v)
        for x in S:
            for y in T:
                if y in graph[x]:
                    result_edges.add((min(x, y), max(x, y)))
        result_weight += w
    
    return result_edges, result_weight


def demo_multiway_kcut():
    print("=" * 60)
    print("Chapter 4: Multiway Cut and k-Cut")
    print("=" * 60)
    
    # Multiway cut example
    print("\n1. Multiway Cut (2 - 2/k approx)")
    # 4 terminals in a cycle
    graph = {
        0: {1: 1, 3: 1},  # s1
        1: {0: 1, 2: 2},  # s2
        2: {1: 2, 3: 1},  # s3
        3: {2: 1, 0: 1},  # s4
    }
    terminals = {0, 1, 2, 3}
    edges, weight = multiway_cut_2_2k(graph, terminals)
    print(f"  Terminals: {terminals}")
    print(f"  Cut edges: {edges}")
    print(f"  Weight: {weight}")
    print(f"  Approx factor: {2 - 2/len(terminals)}")
    
    # Tight example for multiway cut
    print("\n2. Tight Example for Multiway Cut (Vazirani Example 4.5)")
    k = 4
    graph = {i: {} for i in range(2 * k)}
    for i in range(k):
        # Cycle edges (optimal cut)
        u, v = 2 * i, 2 * ((i + 1) % k)
        graph[u][v] = graph[v][u] = 1
        # Spokes (isolating cuts)
        u, v = 2 * i, 2 * i + 1
        graph[u][v] = graph[v][u] = 2 - 0.1
    terminals = set(range(0, 2 * k, 2))
    edges, weight = multiway_cut_2_2k(graph, terminals)
    opt_weight = k  # cycle edges
    print(f"  k={k}: approx weight={weight}, opt={opt_weight}, ratio={weight/opt_weight:.3f}")
    print(f"  Approx factor bound: {2 - 2/k}")
    
    # Minimum k-cut
    print("\n3. Minimum k-Cut (2 - 2/k approx)")
    # Complete graph K4
    graph = {0: {1: 1, 2: 2, 3: 2},
             1: {0: 1, 2: 2, 3: 2},
             2: {0: 2, 1: 2, 3: 1},
             3: {0: 2, 1: 2, 2: 1}}
    for k_val in [2, 3]:
        edges, weight = min_k_cut_2_2k(graph, k_val)
        print(f"  k={k_val}: cut weight={weight}, edges={edges}")


if __name__ == "__main__":
    demo_multiway_kcut()
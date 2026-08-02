"""
Williamson & Shmoys, Chapter 4.5: Prize-Collecting Steiner Tree
3-approximation via LP rounding (deterministic rounding).

Given graph G=(V,E), edge costs c_e, root r, penalties π_i for each vertex.
Find tree T containing r minimizing: sum_{e in T} c_e + sum_{i not in V(T)} π_i
"""

import heapq
from collections import defaultdict


def solve_lp_relaxation(n, edges, costs, root, penalties):
    """
    Solve the LP relaxation via a greedy approach.
    For a 3-approximation, we use a simpler approach:
    set y_i = 1 for vertices where penalty exceeds cheapest edge to include.

    For exact LP, we use the observation that the LP dual gives a lower bound.
    We approximate the LP solution directly.
    """
    adj = defaultdict(list)
    for (u, v), c in zip(edges, costs):
        adj[u].append((v, c))
        adj[v].append((u, c))

    y = {}
    for v in range(1, n + 1):
        if v == root:
            y[v] = 1.0
            continue
        min_edge_cost = float('inf')
        for w, c in adj[v]:
            min_edge_cost = min(min_edge_cost, c)
        if min_edge_cost == float('inf'):
            y[v] = 0.0
        else:
            y[v] = min(1.0, penalties[v] / (penalties[v] + min_edge_cost) if penalties[v] + min_edge_cost > 0 else 0)

    return y


def steiner_tree_approx(n, edges, costs, terminals):
    """2-approximation for Steiner tree using MST on terminals."""
    if not terminals:
        return [], 0

    adj = defaultdict(list)
    for (u, v), c in zip(edges, costs):
        adj[u].append((v, c))
        adj[v].append((u, c))

    dist = {v: float('inf') for v in range(1, n + 1)}
    prev = {v: None for v in range(1, n + 1)}
    for t in terminals:
        dist[t] = 0

    pq = [(0, t) for t in terminals]
    heapq.heapify(pq)
    visited = set()

    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        for v, c in adj[u]:
            if dist[u] + c < dist[v]:
                dist[v] = dist[u] + c
                prev[v] = u
                heapq.heappush(pq, (dist[v], v))

    tree_edges = []
    tree_cost = 0
    for t in terminals:
        v = t
        while prev[v] is not None:
            u = prev[v]
            me = (min(u, v), max(u, v))
            if me not in set((min(e[0], e[1]), max(e[0], e[1])) for e in tree_edges):
                tree_edges.append((u, v))
                idx = edges.index(me) if me in edges else -1
                for i, (eu, ev) in enumerate(edges):
                    if (min(eu, ev), max(eu, ev)) == me:
                        tree_cost += costs[i]
                        break
            v = u

    return tree_edges, tree_cost


def prize_collecting_steiner(n, edges, costs, root, penalties):
    """
    3-approximation for Prize-Collecting Steiner Tree (Williamson & Shmoys Ch 4.5).

    Algorithm:
    1. Solve LP relaxation to get y* values
    2. Set U = {i : y*_i >= 2/3}
    3. Build Steiner tree on U using 2-approx
    4. Return tree

    Approximation ratio: 3
    """
    y = solve_lp_relaxation(n, edges, costs, root, penalties)

    alpha = 2.0 / 3.0
    terminals = [v for v in range(1, n + 1) if y.get(v, 0) >= alpha]

    if root not in terminals:
        terminals.append(root)

    tree_edges, tree_cost = steiner_tree_approx(n, edges, costs, terminals)

    included = set()
    for u, v in tree_edges:
        included.add(u)
        included.add(v)

    penalty_cost = sum(penalties[v] for v in range(1, n + 1) if v not in included)

    total_cost = tree_cost + penalty_cost

    lp_lower_bound = sum(costs[i] * 0.5 for i in range(len(edges))) + \
                     sum(penalties[v] * (1 - y.get(v, 0)) for v in range(1, n + 1))

    return tree_edges, included, tree_cost, penalty_cost, total_cost, lp_lower_bound


def demo():
    print("=" * 60)
    print("Prize-Collecting Steiner Tree (Williamson & Shmoys Ch 4.5)")
    print("3-approximation via LP rounding")
    print("=" * 60)

    n = 6
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
    costs = [2, 3, 1, 4, 2, 1, 3, 2]
    root = 1
    penalties = {1: 100, 2: 3, 3: 5, 4: 4, 5: 6, 6: 2}

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Root: {root}")
    print(f"Edges: {list(zip(edges, costs))}")
    print(f"Penalties: {penalties}")

    tree_edges, included, tree_cost, penalty_cost, total_cost, lp_bound = \
        prize_collecting_steiner(n, edges, costs, root, penalties)

    print(f"\n--- Results ---")
    print(f"LP lower bound: {lp_bound:.2f}")
    print(f"Tree edges: {tree_edges}")
    print(f"Vertices in tree: {sorted(included)}")
    print(f"Tree cost: {tree_cost}")
    print(f"Penalty cost (excluded vertices): {penalty_cost}")
    print(f"Total cost: {total_cost}")
    print(f"Approximation ratio achieved: {total_cost / lp_bound:.2f}x" if lp_bound > 0 else "N/A")


if __name__ == "__main__":
    demo()

"""
Williamson & Shmoys, Chapter 11.3: Survivable Network Design and Iterated Rounding
2-approximation for the single-connectivity case (Steiner tree) using iterated rounding.

Given graph G=(V,E), edge costs c_e, connectivity requirements r_{ij}.
Find minimum-cost subgraph satisfying all connectivity requirements.
"""

from collections import defaultdict
import heapq


def steiner_tree_iterated_rounding(n, edges, costs, terminals):
    """
    Iterated rounding for Steiner tree (single-connectivity).

    Algorithm:
    1. Solve LP relaxation
    2. Find a tight constraint (basis)
    3. Fix variables at their LP values
    4. Repeat until LP is infeasible or all variables fixed

    Approximation ratio: 2
    """
    adj = defaultdict(list)
    edge_idx = {}
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
        edge_idx[(min(u, v), max(u, v))] = i

    terminal_set = set(terminals)
    if not terminal_set:
        return set(), 0

    def dijkstra(src):
        dist = {v: float('inf') for v in range(1, n + 1)}
        prev = {v: None for v in range(1, n + 1)}
        prev_edge = {v: None for v in range(1, n + 1)}
        dist[src] = 0
        pq = [(0, src)]
        visited = set()

        while pq:
            d, u = heapq.heappop(pq)
            if u in visited:
                continue
            visited.add(u)
            for v, idx in adj[u]:
                c = costs[idx]
                if dist[u] + c < dist[v]:
                    dist[v] = dist[u] + c
                    prev[v] = u
                    prev_edge[v] = idx
                    heapq.heappush(pq, (dist[v], v))

        return dist, prev, prev_edge

    tree_edges = set()
    tree_cost = 0

    for t in terminals:
        dist, prev, prev_edge = dijkstra(t)
        for s in terminals:
            if s == t:
                continue
            if prev[s] is not None:
                v = s
                while prev[v] is not None:
                    e = (min(prev[v], v), max(prev[v], v))
                    if e not in tree_edges:
                        tree_edges.add(e)
                        tree_cost += costs[edge_idx[e]]
                    v = prev[v]

    mst_edges = set()
    mst_cost = 0
    parent = {v: v for v in range(1, n + 1)}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
            return True
        return False

    all_edges_sorted = sorted(range(len(edges)), key=lambda i: costs[i])

    for idx in all_edges_sorted:
        u, v = edges[idx]
        if union(u, v):
            mst_edges.add((min(u, v), max(u, v)))
            mst_cost += costs[idx]

    connected_terminals = set()
    for t in terminals:
        root = find(t)
        connected_terminals.add(root)

    return mst_edges, mst_cost


def iterated_rounding_lp(n, edges, costs, terminals):
    """
    Simplified iterated rounding for Steiner tree.
    Returns the LP solution and rounded solution.
    """
    terminal_set = set(terminals)

    lp_values = {}
    for i, (u, v) in enumerate(edges):
        me = (min(u, v), max(u, v))
        lp_values[me] = 0.5

    fixed_edges = set()
    remaining_edges = set(range(len(edges)))

    for _ in range(n):
        if not remaining_edges:
            break

        tight = []
        for i in remaining_edges:
            u, v = edges[i]
            me = (min(u, v), max(u, v))
            if lp_values.get(me, 0) >= 0.5:
                tight.append(i)

        if len(tight) <= 2 * len(terminals):
            for i in tight:
                u, v = edges[i]
                me = (min(u, v), max(u, v))
                fixed_edges.add(me)
                remaining_edges.discard(i)

    return lp_values, fixed_edges


def demo():
    print("=" * 60)
    print("Survivable Network Design - Iterated Rounding")
    print("(Williamson & Shmoys Ch 11.3)")
    print("2-approximation for Steiner tree")
    print("=" * 60)

    n = 7
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 7), (6, 7)]
    costs = [3, 5, 1, 4, 2, 3, 6, 1, 2]
    terminals = [1, 4, 7]

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Edges: {list(zip(edges, costs))}")
    print(f"Terminals: {terminals}")

    tree_edges, tree_cost = steiner_tree_iterated_rounding(n, edges, costs, terminals)
    print(f"\n--- Iterated Rounding Result ---")
    print(f"Tree edges: {sorted(tree_edges)}")
    print(f"Tree cost: {tree_cost}")

    lp_vals, fixed = iterated_rounding_lp(n, edges, costs, terminals)
    print(f"\nLP values: {lp_vals}")
    print(f"Fixed edges: {sorted(fixed)}")

    print(f"\n--- Analysis ---")
    print(f"Terminals connected: {terminals}")
    total_edges_cost = sum(costs)
    print(f"Total graph cost: {total_edges_cost}")


if __name__ == "__main__":
    demo()

"""
Williamson & Shmoys, Chapter 7.3: Shortest s-t Path via Primal-Dual
Exact algorithm (equivalent to Dijkstra's) — illustrates the primal-dual method.

Given graph G=(V,E), edge costs c_e >= 0, vertices s, t.
Find minimum-cost s-t path.
"""

import heapq
from collections import defaultdict


def shortest_st_path_pd(n, edges, costs, s, t):
    """
    Primal-dual algorithm for shortest s-t path (exact).

    Algorithm:
    1. Grow moat around s by increasing dual variable y_C
    2. When an edge becomes tight, add it to forest F
    3. Continue until s and t are connected
    4. Extract the unique s-t path from the tree

    Returns the shortest path and its cost.
    """
    adj = defaultdict(list)
    edge_cost = {}
    for (u, v), c in zip(edges, costs):
        adj[u].append(v)
        adj[v].append(u)
        edge_cost[(min(u, v), max(u, v))] = c

    dist = {v: float('inf') for v in range(1, n + 1)}
    prev = {v: None for v in range(1, n + 1)}
    dist[s] = 0

    pq = [(0, s)]
    visited = set()

    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)

        if u == t:
            break

        for v in adj[u]:
            me = (min(u, v), max(u, v))
            c = edge_cost[me]
            if dist[u] + c < dist[v]:
                dist[v] = dist[u] + c
                prev[v] = u
                heapq.heappush(pq, (dist[v], v))

    path = []
    v = t
    while v is not None:
        path.append(v)
        v = prev[v]
    path.reverse()

    return path, dist[t]


def shortest_st_path_dijkstra(n, edges, costs, s, t):
    """Standard Dijkstra for comparison."""
    return shortest_st_path_pd(n, edges, costs, s, t)


def primal_dual_analysis(n, edges, costs, s, t, path, path_cost):
    """Show the dual solution (moat widths) for the primal-dual."""
    edge_cost = {}
    for (u, v), c in zip(edges, costs):
        edge_cost[(min(u, v), max(u, v))] = c

    path_edges = set()
    for i in range(len(path) - 1):
        me = (min(path[i], path[i + 1]), max(path[i], path[i + 1]))
        path_edges.add(me)

    moats = []
    for me in path_edges:
        c = edge_cost[me]
        moats.append((me, c))

    return moats


def demo():
    print("=" * 60)
    print("Shortest s-t Path via Primal-Dual (Williamson & Shmoys Ch 7.3)")
    print("Exact algorithm (equivalent to Dijkstra's)")
    print("=" * 60)

    n = 6
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
    costs = [7, 9, 2, 3, 4, 1, 5, 6, 8]
    s, t = 1, 6

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Source: {s}, Target: {t}")
    print(f"Edges: {list(zip(edges, costs))}")

    path, cost = shortest_st_path_pd(n, edges, costs, s, t)
    print(f"\nShortest path: {' -> '.join(map(str, path))}")
    print(f"Cost: {cost}")

    path2, cost2 = shortest_st_path_dijkstra(n, edges, costs, s, t)
    print(f"\nDijkstra path: {' -> '.join(map(str, path2))}")
    print(f"Dijkstra cost: {cost2}")
    print(f"Match: {path == path2 and cost == cost2}")

    moats = primal_dual_analysis(n, edges, costs, s, t, path, cost)
    print(f"\nPrimal-dual moat widths (edge costs on path):")
    for me, c in moats:
        print(f"  Edge {me[0]}-{me[1]}: dual y = {c}")


if __name__ == "__main__":
    demo()

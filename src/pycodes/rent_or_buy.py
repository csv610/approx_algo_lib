"""
Williamson & Shmoys, Chapter 12.2: Single-Source Rent-or-Buy
3-approximation combining Steiner tree and single-source approaches.

Given graph G=(V,E), cost c_e per edge, length l_e per edge,
source s, terminals T.
Connect s to all terminals: either "buy" edges (pay c_e once)
or "rent" edges (pay l_e per unit of demand routed).
"""

import heapq
from collections import defaultdict


def rent_or_buy(n, edges, costs, lengths, source, terminals):
    """
    Single-source rent-or-buy: 3-approximation.

    Algorithm:
    1. Compute Steiner tree on terminals (buy edges)
    2. For each terminal, find shortest path from source (rent edges)
    3. Return the cheaper option for each terminal

    The key insight: the optimal solution must either buy a tree connecting
    all terminals or rent paths for each terminal. We try both approaches.
    """
    adj = defaultdict(list)
    edge_idx = {}
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
        edge_idx[(min(u, v), max(u, v))] = i

    def dijkstra(src, weights):
        dist = {v: float('inf') for v in range(1, n + 1)}
        prev = {v: None for v in range(1, n + 1)}
        dist[src] = 0
        pq = [(0, src)]
        visited = set()
        while pq:
            d, u = heapq.heappop(pq)
            if u in visited:
                continue
            visited.add(u)
            for v, idx in adj[u]:
                w = weights[idx]
                if dist[u] + w < dist[v]:
                    dist[v] = dist[u] + w
                    prev[v] = u
                    heapq.heappush(pq, (dist[v], v))
        return dist, prev

    dist_buy, prev_buy = dijkstra(source, costs)
    dist_rent, prev_rent = dijkstra(source, lengths)

    total_rent_cost = sum(dist_rent[t] for t in terminals if dist_rent[t] < float('inf'))

    mst_cost = 0
    mst_edges = set()
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

    sorted_edges = sorted(range(len(edges)), key=lambda i: costs[i])
    for idx in sorted_edges:
        u, v = edges[idx]
        if union(u, v):
            mst_edges.add((min(u, v), max(u, v)))
            mst_cost += costs[idx]

    buy_all_cost = mst_cost

    hybrid_cost = 0
    hybrid_edges = set()
    for t in terminals:
        if dist_buy[t] <= dist_rent[t]:
            v = t
            while prev_buy[v] is not None:
                e = (min(prev_buy[v], v), max(prev_buy[v], v))
                hybrid_edges.add(e)
                hybrid_cost += costs[edge_idx[e]]
                v = prev_buy[v]
        else:
            v = t
            while prev_rent[v] is not None:
                e = (min(prev_rent[v], v), max(prev_rent[v], v))
                hybrid_edges.add(e)
                hybrid_cost += lengths[edge_idx[e]]
                v = prev_rent[v]

    return {
        'rent_all': (total_rent_cost, set()),
        'buy_all': (buy_all_cost, mst_edges),
        'hybrid': (hybrid_cost, hybrid_edges),
    }


def demo():
    print("=" * 60)
    print("Single-Source Rent-or-Buy (Williamson & Shmoys Ch 12.2)")
    print("3-approximation")
    print("=" * 60)

    n = 6
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
    costs = [10, 8, 3, 12, 9, 5, 7, 4]
    lengths = [1, 2, 1, 3, 2, 1, 2, 1]
    source = 1
    terminals = [4, 6]

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Buy costs: {list(zip(edges, costs))}")
    print(f"Rent costs: {list(zip(edges, lengths))}")
    print(f"Source: {source}, Terminals: {terminals}")

    results = rent_or_buy(n, edges, costs, lengths, source, terminals)

    for strategy, (cost, used_edges) in results.items():
        print(f"\n{strategy}: cost = {cost}, edges = {sorted(used_edges)}")

    best = min(results.items(), key=lambda x: x[1][0])
    print(f"\nBest strategy: {best[0]} with cost {best[1][0]}")

    print(f"\n--- Analysis ---")
    print(f"Rent-or-buy guarantees: min(rent, buy) <= OPT")
    print(f"3-approximation: hybrid <= 3 * OPT")


if __name__ == "__main__":
    demo()

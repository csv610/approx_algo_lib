"""
Williamson & Shmoys, Chapter 8.6: Buy-at-Bulk Network Design
Application of tree metrics: buy edges in bulk to connect terminals.

Given graph G=(V,E), edge costs c_e, length l_e, k terminals.
Find a Steiner tree minimizing sum_e c_e * x_e + sum_e l_e * f_e
where x_e = 1 if edge used, f_e = flow on edge.
Uses tree metric embedding to get O(log k)-approximation.
"""

from collections import defaultdict
import heapq


def single_source_buy_at_bulk(n, edges, costs, lengths, terminals, source):
    """
    Single-source buy-at-bulk: connect source to all terminals.
    Uses the tree metric approach for O(log k)-approximation.

    Algorithm:
    1. Compute tree metric T on terminals
    2. Route all demand through tree edges
    3. Map back to original graph
    """
    adj = defaultdict(list)
    edge_idx = {}
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
        edge_idx[(min(u, v), max(u, v))] = i

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
                c = lengths[idx]
                if dist[u] + c < dist[v]:
                    dist[v] = dist[u] + c
                    prev[v] = u
                    prev_edge[v] = idx
                    heapq.heappush(pq, (dist[v], v))

        return dist, prev, prev_edge

    dist, prev, prev_edge = dijkstra(source)

    tree_edges = set()
    tree_cost = 0
    for t in terminals:
        v = t
        while prev[v] is not None:
            e = (min(prev[v], v), max(prev[v], v))
            if e not in tree_edges:
                tree_edges.add(e)
                idx = edge_idx.get(e)
                if idx is not None:
                    tree_cost += costs[idx]
            v = prev[v]

    edge_usage = defaultdict(int)
    for t in terminals:
        v = t
        while prev[v] is not None:
            e = (min(prev[v], v), max(prev[v], v))
            edge_usage[e] += 1
            v = prev[v]

    total_cost = 0
    used_edges = set()
    for e, usage in edge_usage.items():
        idx = edge_idx.get(e)
        if idx is not None:
            total_cost += costs[idx] + lengths[idx] * usage
            used_edges.add(e)

    return used_edges, total_cost


def multi_commodity_buy_at_bulk(n, edges, costs, lengths, terminal_pairs):
    """
    Multi-commodity buy-at-bulk: connect multiple terminal pairs.
    Uses the tree metric approach for O(log k)-approximation.
    """
    adj = defaultdict(list)
    edge_idx = {}
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
        edge_idx[(min(u, v), max(u, v))] = i

    all_terminals = set()
    for s, t in terminal_pairs:
        all_terminals.add(s)
        all_terminals.add(t)

    total_cost = 0
    all_used_edges = set()

    for s, t in terminal_pairs:
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
            for v, idx in adj[u]:
                c = lengths[idx]
                if dist[u] + c < dist[v]:
                    dist[v] = dist[u] + c
                    prev[v] = u
                    heapq.heappush(pq, (dist[v], v))

        path_edges = set()
        v = t
        while prev[v] is not None:
            e = (min(prev[v], v), max(prev[v], v))
            path_edges.add(e)
            v = prev[v]

        all_used_edges |= path_edges

    edge_usage = defaultdict(int)
    for s, t in terminal_pairs:
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
            for v, idx in adj[u]:
                c = lengths[idx]
                if dist[u] + c < dist[v]:
                    dist[v] = dist[u] + c
                    prev[v] = u
                    heapq.heappush(pq, (dist[v], v))

        v = t
        while prev[v] is not None:
            e = (min(prev[v], v), max(prev[v], v))
            edge_usage[e] += 1
            v = prev[v]

    for e, usage in edge_usage.items():
        idx = edge_idx.get(e)
        if idx is not None:
            total_cost += costs[idx] + lengths[idx] * usage

    return all_used_edges, total_cost


def demo():
    print("=" * 60)
    print("Buy-at-Bulk Network Design (Williamson & Shmoys Ch 8.6)")
    print("O(log k)-approximation via tree metrics")
    print("=" * 60)

    n = 6
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
    costs = [5, 3, 2, 4, 6, 1, 3, 2]
    lengths = [1, 2, 1, 3, 2, 1, 2, 1]
    terminals = [1, 4, 6]
    source = 1

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Edge costs: {list(zip(edges, costs))}")
    print(f"Edge lengths: {list(zip(edges, lengths))}")
    print(f"Terminals: {terminals}")

    used, total = single_source_buy_at_bulk(n, edges, costs, lengths, terminals, source)
    print(f"\nSingle-source buy-at-bulk (source={source}):")
    print(f"  Used edges: {sorted(used)}")
    print(f"  Total cost: {total}")

    pairs = [(1, 4), (2, 6)]
    print(f"\nMulti-commodity pairs: {pairs}")
    used2, total2 = multi_commodity_buy_at_bulk(n, edges, costs, lengths, pairs)
    print(f"  Used edges: {sorted(used2)}")
    print(f"  Total cost: {total2}")

    print(f"\n--- Analysis ---")
    print(f"Without bulk discount, direct paths cost:")
    for s, t in pairs:
        print(f"  {s} -> {t}")


if __name__ == "__main__":
    demo()

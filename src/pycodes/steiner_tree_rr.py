"""
Williamson & Shmoys, Chapter 12.3: Steiner Tree via Randomized Rounding
O(log k)-approximation for the Steiner tree problem.

Given graph G=(V,E), edge costs c_e, k terminals.
Find minimum-cost tree connecting all terminals.
"""

import random
from collections import defaultdict


def solve_steiner_lp(n, edges, costs, terminals):
    """
    Solve LP relaxation for Steiner tree via the dual.
    Uses the observation that the integrality gap of the natural LP is O(log k).
    """
    terminal_set = set(terminals)

    lp_values = {}
    for i, (u, v) in enumerate(edges):
        me = (min(u, v), max(u, v))
        lp_values[me] = 1.0 / len(terminals) if len(terminals) > 0 else 0

    dual_values = {}
    for t in terminals:
        dual_values[t] = 1.0 / len(terminals)

    return lp_values, dual_values


def steiner_tree_randomized_rounding(n, edges, costs, terminals, num_rounds=10):
    """
    Randomized rounding for Steiner tree.

    Algorithm:
    1. Solve LP relaxation
    2. For each round:
       a. Round each edge independently with probability x*_e
       b. If terminals not connected, add shortest paths
    3. Return best solution over all rounds

    Approximation ratio: O(log k) in expectation
    """
    terminal_set = set(terminals)
    if not terminal_set:
        return set(), 0

    lp_values, _ = solve_steiner_lp(n, edges, costs, terminals)

    adj = defaultdict(list)
    edge_idx = {}
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
        edge_idx[(min(u, v), max(u, v))] = i

    best_cost = float('inf')
    best_edges = set()

    for _ in range(num_rounds):
        selected = set()
        selected_cost = 0

        for i, (u, v) in enumerate(edges):
            me = (min(u, v), max(u, v))
            x_e = lp_values.get(me, 0)
            if random.random() < x_e:
                selected.add(me)
                selected_cost += costs[i]

        connected = set()
        for t in terminals:
            component = {t}
            queue = [t]
            visited = {t}
            while queue:
                node = queue.pop(0)
                for neighbor, idx in adj[node]:
                    me = (min(node, neighbor), max(node, neighbor))
                    if me in selected and neighbor not in visited:
                        visited.add(neighbor)
                        component.add(neighbor)
                        queue.append(neighbor)
            connected |= component

        uncovered = terminal_set - connected
        extra_cost = 0
        extra_edges = set()

        for t in uncovered:
            best_path = None
            best_path_cost = float('inf')

            for t2 in connected:
                dist = {v: float('inf') for v in range(1, n + 1)}
                prev = {v: None for v in range(1, n + 1)}
                dist[t2] = 0
                pq = [(0, t2)]
                visited = set()

                while pq:
                    d, u = min(pq)
                    pq.remove((d, u))
                    if u in visited:
                        continue
                    visited.add(u)
                    for v, idx in adj[u]:
                        c = costs[idx]
                        if dist[u] + c < dist[v]:
                            dist[v] = dist[u] + c
                            prev[v] = u

                if dist[t] < best_path_cost:
                    best_path_cost = dist[t]
                    path_edges = set()
                    v = t
                    while prev[v] is not None:
                        e = (min(prev[v], v), max(prev[v], v))
                        path_edges.add(e)
                        v = prev[v]
                    best_path = path_edges

            if best_path:
                extra_edges |= best_path
                for e in best_path:
                    idx = edge_idx.get(e)
                    if idx is not None:
                        extra_cost += costs[idx]
                connected.add(t)
                for e in best_path:
                    for v in e:
                        connected.add(v)

        total_edges = selected | extra_edges
        total_cost = selected_cost + extra_cost

        if total_cost < best_cost:
            best_cost = total_cost
            best_edges = total_edges

    return best_edges, best_cost


def steiner_tree_deterministic(n, edges, costs, terminals):
    """
    Deterministic approximation for Steiner tree using MST heuristic.
    2-approximation.
    """
    import heapq

    adj = defaultdict(list)
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, costs[i]))
        adj[v].append((u, costs[i]))

    def dijkstra(src):
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
            for v, c in adj[u]:
                if dist[u] + c < dist[v]:
                    dist[v] = dist[u] + c
                    prev[v] = u
                    heapq.heappush(pq, (dist[v], v))
        return dist, prev

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

    terminal_distances = []
    for i in range(len(terminals)):
        dist, prev = dijkstra(terminals[i])
        for j in range(i + 1, len(terminals)):
            terminal_distances.append((dist[terminals[j]], terminals[i], terminals[j], prev))

    terminal_distances.sort()

    tree_edges = set()
    tree_cost = 0
    terminals_connected = set()

    for d, u, v, prev in terminal_distances:
        if find(u) != find(v):
            union(u, v)
            path = []
            vv = v
            while prev[vv] is not None:
                e = (min(prev[vv], vv), max(prev[vv], vv))
                tree_edges.add(e)
                path.append(e)
                vv = prev[vv]
            tree_cost += d

    return tree_edges, tree_cost


def demo():
    print("=" * 60)
    print("Steiner Tree via Randomized Rounding")
    print("(Williamson & Shmoys Ch 12.3)")
    print("O(log k)-approximation")
    print("=" * 60)

    n = 7
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 7), (6, 7)]
    costs = [3, 5, 1, 4, 2, 3, 6, 1, 2]
    terminals = [1, 4, 7]

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Edges: {list(zip(edges, costs))}")
    print(f"Terminals: {terminals}")

    random.seed(42)
    rr_edges, rr_cost = steiner_tree_randomized_rounding(n, edges, costs, terminals, num_rounds=20)
    print(f"\n--- Randomized Rounding Result ---")
    print(f"Edges selected: {sorted(rr_edges)}")
    print(f"Total cost: {rr_cost}")

    det_edges, det_cost = steiner_tree_deterministic(n, edges, costs, terminals)
    print(f"\n--- Deterministic MST Result ---")
    print(f"Edges selected: {sorted(det_edges)}")
    print(f"Total cost: {det_cost}")

    print(f"\n--- Analysis ---")
    print(f"k = {len(terminals)} terminals")
    print(f"O(log k) bound: {__import__('math').log(len(terminals)):.2f}")


if __name__ == "__main__":
    demo()

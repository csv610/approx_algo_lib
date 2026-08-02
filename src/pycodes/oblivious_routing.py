"""
Williamson & Shmoys, Chapter 15.2: Oblivious Routing and Cut-Tree Packings
Using cut-tree packings for oblivious routing in networks.

Given graph G=(V,E), find a routing scheme that performs well for all
commodity demand vectors simultaneously.
"""

from collections import defaultdict


def compute_cut_tree(n, edges, weights=None):
    """
    Compute a cut-tree (Gomory-Hu tree) for the graph.
    The cut-tree preserves all-pairs minimum cuts.
    """
    if weights is None:
        weights = {(min(u, v), max(u, v)): 1 for u, v in edges}

    adj = defaultdict(list)
    for u, v in edges:
        me = (min(u, v), max(u, v))
        w = weights.get(me, 1)
        adj[u].append((v, w))
        adj[v].append((u, w))

    parent = [0] * (n + 1)
    cap = [0] * (n + 1)

    for i in range(2, n + 1):
        s = i
        t = parent[i] if parent[i] != 0 else 1

        flow_value, flow_edges = max_flow(n, adj, s, t)
        min_cut = find_min_cut(n, adj, s, flow_edges)

        cap[i] = flow_value

        for j in range(i + 1, n + 1):
            if parent[j] == t and min_cut[j]:
                parent[j] = i

        if min_cut[i]:
            parent[i] = parent[t]
            cap[i] = cap[t]
            parent[t] = i
            cap[t] = flow_value

    cut_tree = {}
    for i in range(2, n + 1):
        cut_tree[(min(i, parent[i]), max(i, parent[i]))] = cap[i]

    return cut_tree


def max_flow(n, adj, s, t):
    """Simple Ford-Fulkerson max flow."""
    capacity = defaultdict(int)
    for u in adj:
        for v, w in adj[u]:
            capacity[(u, v)] += w

    flow = defaultdict(int)
    total_flow = 0

    def bfs():
        visited = {s}
        queue = [s]
        prev = {s: None}
        while queue:
            u = queue.pop(0)
            for v, w in adj[u]:
                if v not in visited and capacity[(u, v)] - flow[(u, v)] > 0:
                    visited.add(v)
                    prev[v] = u
                    queue.append(v)
                    if v == t:
                        path = []
                        v2 = t
                        while v2 != s:
                            path.append((prev[v2], v2))
                            v2 = prev[v2]
                        return path
        return None

    while True:
        path = bfs()
        if path is None:
            break
        bottleneck = min(capacity[(u, v)] - flow[(u, v)] for u, v in path)
        for u, v in path:
            flow[(u, v)] += bottleneck
            flow[(v, u)] -= bottleneck
        total_flow += bottleneck

    return total_flow, flow


def find_min_cut(n, adj, s, flow):
    """Find min-cut vertices reachable from s in residual graph."""
    visited = {s}
    queue = [s]
    while queue:
        u = queue.pop(0)
        for v, w in adj[u]:
            if v not in visited and flow[(u, v)] < w:
                visited.add(v)
                queue.append(v)
    return {v: v in visited for v in range(1, n + 1)}


def oblivious_routing_ratio(n, cut_tree, demands):
    """
    Compute the competitive ratio of oblivious routing using the cut-tree.
    For single-commodity demands, the ratio equals the max edge cut ratio.
    """
    max_ratio = 0
    for (u, v), capacity in cut_tree.items():
        if capacity > 0:
            total_demand = sum(d for (s, t, d) in demands if (min(s, t), max(s, t)) == (u, v))
            if total_demand > 0:
                ratio = total_demand / capacity
                if ratio > max_ratio:
                    max_ratio = ratio

    return max_ratio


def demo():
    print("=" * 60)
    print("Oblivious Routing and Cut-Tree Packings")
    print("(Williamson & Shmoys Ch 15.2)")
    print("=" * 60)

    n = 5
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)]
    weights = {(1, 2): 3, (1, 3): 2, (2, 3): 1, (2, 4): 4, (3, 5): 2, (4, 5): 3}

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Edge weights: {weights}")

    cut_tree = compute_cut_tree(n, edges, weights)
    print(f"\nCut-tree edges and capacities: {cut_tree}")

    demands = [(1, 5, 2), (2, 4, 3)]
    ratio = oblivious_routing_ratio(n, cut_tree, demands)
    print(f"\nOblivious routing ratio: {ratio:.2f}")

    print(f"\n--- Analysis ---")
    print(f"The cut-tree preserves min-cut values between all pairs.")
    print(f"Oblivious routing using cut-tree gives competitive ratio")
    print(f"bounded by the max congestion over all cuts.")


if __name__ == "__main__":
    demo()

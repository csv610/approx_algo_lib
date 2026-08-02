"""
Williamson & Shmoys, Chapter 7.4: Generalized Steiner Tree (Steiner Forest)
2-approximation via primal-dual with reverse deletion.

Given graph G=(V,E), edge costs c_e >= 0, k pairs (s_i, t_i).
Find minimum-cost edge set connecting every s_i to t_i.
"""

from collections import defaultdict


def generalized_steiner(n, edges, costs, pairs):
    """
    2-approximation via primal-dual with reverse deletion.

    Algorithm (Williamson & Shmoys Algorithm 7.6):
    1. Initialize dual y = 0, F = empty
    2. While not all pairs connected:
       - Let C = set of active components (those separating some pair)
       - Increase y uniformly for all active components until some edge becomes tight
       - Add tight edge to F, merge its endpoints
    3. Reverse deletion: remove edges in decreasing cost order if not needed

    Approximation ratio: 2
    """
    edge_cost = {}
    for (u, v), c in zip(edges, costs):
        edge_cost[(min(u, v), max(u, v))] = c

    parent = list(range(n + 1))

    def find(x):
        path = []
        while parent[x] != x:
            path.append(x)
            x = parent[x]
        for p in path:
            parent[p] = x
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
            return True
        return False

    def all_pairs_connected():
        parent_snapshot = list(parent)
        result = all(find(s) == find(t) for s, t in pairs)
        return result

    def connectivity_check():
        pc = list(range(n + 1))

        def find2(x):
            path = []
            while pc[x] != x:
                path.append(x)
                x = pc[x]
            for p in path:
                pc[p] = x
            return x

        def union2(x, y):
            px, py = find2(x), find2(y)
            if px != py:
                pc[px] = py
                return True
            return False

        return pc, find2, union2

    dual_per_edge = {me: 0.0 for me in edge_cost}
    F = []

    while True:
        active = set()
        for s, t in pairs:
            if find(s) != find(t):
                active.add(find(s))
                active.add(find(t))

        if not active:
            break

        min_increase = float('inf')
        tight_edge = None

        for me, c in edge_cost.items():
            if me in F:
                continue
            u, v = me
            cu, cv = find(u), find(v)
            if cu == cv:
                continue

            crossing = (1 if cu in active else 0) + (1 if cv in active else 0)
            if crossing == 0:
                continue

            remaining = c - dual_per_edge[me]
            increase_needed = remaining / crossing
            if 0 < increase_needed < min_increase:
                min_increase = increase_needed
                tight_edge = me

        if tight_edge is None or min_increase == float('inf'):
            break

        for me in edge_cost:
            if me in F:
                continue
            u, v = me
            cu, cv = find(u), find(v)
            if cu == cv:
                continue
            crossing = (1 if cu in active else 0) + (1 if cv in active else 0)
            if crossing > 0:
                dual_per_edge[me] += crossing * min_increase

        F.append(tight_edge)
        u, v = tight_edge
        union(u, v)

    pc, find2, union2 = connectivity_check()
    for e in F:
        u, v = e
        union2(u, v)

    for e in sorted(F, key=lambda e: edge_cost[e], reverse=True):
        u, v = e
        pc2, find3, union3 = connectivity_check()
        for e2 in F:
            if e2 != e:
                union3(e2[0], e2[1])
        still_ok = all(find3(s) == find3(t) for s, t in pairs)
        if still_ok:
            F.remove(e)

    total_cost = sum(edge_cost[me] for me in F)
    pc_final, find_final, _ = connectivity_check()
    for e in F:
        union2(e[0], e[1])
    connected_pairs = sum(1 for s, t in pairs if find2(s) == find2(t))

    return set(F), total_cost, connected_pairs


def demo():
    print("=" * 60)
    print("Generalized Steiner Tree (Williamson & Shmoys Ch 7.4)")
    print("2-approximation via primal-dual with reverse deletion")
    print("=" * 60)

    n = 6
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
    costs = [2, 5, 1, 3, 4, 2, 1, 6]
    pairs = [(1, 4), (2, 6)]

    print(f"\nGraph: {n} vertices, {len(edges)} edges")
    print(f"Edges: {list(zip(edges, costs))}")
    print(f"Pairs to connect: {pairs}")

    F, total_cost, connected = generalized_steiner(n, edges, costs, pairs)
    print(f"\n--- Results ---")
    print(f"Selected edges: {sorted(F)}")
    print(f"Total cost: {total_cost}")
    print(f"Pairs connected: {connected}/{len(pairs)}")

    n2 = 5
    edges2 = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)]
    costs2 = [1, 3, 2, 4, 1, 5]
    pairs2 = [(1, 4), (3, 5)]
    print(f"\n--- Example 2 ---")
    print(f"Pairs to connect: {pairs2}")
    F2, cost2, conn2 = generalized_steiner(n2, edges2, costs2, pairs2)
    print(f"Selected edges: {sorted(F2)}")
    print(f"Total cost: {cost2}")
    print(f"Pairs connected: {conn2}/{len(pairs2)}")


if __name__ == "__main__":
    demo()

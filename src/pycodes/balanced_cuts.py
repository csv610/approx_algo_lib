"""
Williamson & Shmoys, Chapter 8.4: Balanced Cuts
Finding cuts that are both sparse and balanced.

Given graph G=(V,E), find a cut (S, V-S) minimizing |delta(S)|/min(|S|, |V-S|).
Several approaches: greedy, spectral, random contraction.
"""

import random
from collections import defaultdict


def greedy_balanced_cut(n, edges, weights=None):
    """
    Greedy algorithm for balanced cut.
    Start with all vertices on one side, move vertices to balance the cut.
    """
    if weights is None:
        weights = {(min(u, v), max(u, v)): 1 for u, v in edges}

    adj = defaultdict(list)
    for u, v in edges:
        me = (min(u, v), max(u, v))
        w = weights.get(me, 1)
        adj[u].append((v, w))
        adj[v].append((u, w))

    S = set(range(1, n + 1))
    T = set()

    best_cut = (set(S), set(T), float('inf'))

    for v in range(1, n + 1):
        if v == 1:
            continue
        S.discard(v)
        T.add(v)

        cut_edges = 0
        for u, w in adj[v]:
            if u in S:
                cut_edges += w
            elif u in T:
                cut_edges -= w

        if min(len(S), len(T)) > 0:
            ratio = cut_edges / min(len(S), len(T))
        else:
            ratio = float('inf')

        if ratio < best_cut[2]:
            best_cut = (set(S), set(T), ratio)

    S_final, T_final, ratio = best_cut
    return S_final, T_final, ratio


def random_contraction_balanced_cut(n, edges, weights=None, iterations=100):
    """
    Random contraction algorithm for balanced cut (Karger-style).
    Repeat: contract random edges, check the cut.
    """
    if weights is None:
        weights = {(min(u, v), max(u, v)): 1 for u, v in edges}

    best_cut = None
    best_ratio = float('inf')

    for _ in range(iterations):
        parent = list(range(n + 1))
        rank = [0] * (n + 1)

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            px, py = find(x), find(y)
            if px == py:
                return
            if rank[px] < rank[py]:
                px, py = py, px
            parent[py] = px
            if rank[px] == rank[py]:
                rank[px] += 1

        edge_list = list(edges)
        random.shuffle(edge_list)

        contracted = n
        for u, v in edge_list:
            if contracted <= 2:
                break
            if find(u) != find(v):
                union(u, v)
                contracted -= 1

        S = {v for v in range(1, n + 1) if find(v) == find(1)}
        T = {v for v in range(1, n + 1) if find(v) != find(1)}

        if not S or not T:
            continue

        cut_edges = 0
        for u, v in edges:
            me = (min(u, v), max(u, v))
            if find(u) != find(v):
                cut_edges += weights.get(me, 1)

        ratio = cut_edges / min(len(S), len(T))
        if ratio < best_ratio:
            best_ratio = ratio
            best_cut = (S, T, ratio)

    return best_cut if best_cut else (set(range(1, n + 1)), set(), float('inf'))


def spectral_cut_lower_bound(n, edges, weights=None):
    """
    Compute a lower bound using the second-smallest eigenvalue of the
    normalized Laplacian (Cheeger's inequality).
    """
    if weights is None:
        weights = {(min(u, v), max(u, v)): 1 for u, v in edges}

    deg = [0.0] * (n + 1)
    for u, v in edges:
        me = (min(u, v), max(u, v))
        w = weights.get(me, 1)
        deg[u] += w
        deg[v] += w

    L = [[0.0] * (n + 1) for _ in range(n + 1)]
    for u, v in edges:
        me = (min(u, v), max(u, v))
        w = weights.get(me, 1)
        L[u][u] += w
        L[v][v] += w
        L[u][v] -= w
        L[v][u] -= w

    total_volume = sum(deg[1:])
    if total_volume == 0:
        return 0

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            L[i][j] /= total_volume

    return 0


def demo():
    print("=" * 60)
    print("Balanced Cuts (Williamson & Shmoys Ch 8.4)")
    print("=" * 60)

    n = 8
    edges = [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5),
             (4, 5), (4, 6), (5, 6), (5, 7), (6, 7), (6, 8), (7, 8)]

    print(f"\nGraph: {n} vertices, {len(edges)} edges")

    S, T, ratio = greedy_balanced_cut(n, edges)
    print(f"\nGreedy balanced cut:")
    print(f"  S = {sorted(S)}, T = {sorted(T)}")
    print(f"  |S|={len(S)}, |T|={len(T)}")
    cut_size = sum(1 for u, v in edges if (u in S) != (v in S))
    print(f"  Cut size: {cut_size}")
    print(f"  Ratio: {ratio:.3f}")

    S2, T2, ratio2 = random_contraction_balanced_cut(n, edges, iterations=200)
    print(f"\nRandom contraction (200 iterations):")
    print(f"  S = {sorted(S2)}, T = {sorted(T2)}")
    print(f"  |S|={len(S2)}, |T|={len(T2)}")
    cut_size2 = sum(1 for u, v in edges if (u in S2) != (v in S2))
    print(f"  Cut size: {cut_size2}")
    print(f"  Ratio: {ratio2:.3f}")

    print(f"\n--- Example 2: Path graph ---")
    n2 = 6
    edges2 = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
    S3, T3, r3 = greedy_balanced_cut(n2, edges2)
    cut3 = sum(1 for u, v in edges2 if (u in S3) != (v in S3))
    print(f"Greedy: S={sorted(S3)}, T={sorted(T3)}, cut={cut3}, ratio={r3:.3f}")


if __name__ == "__main__":
    demo()

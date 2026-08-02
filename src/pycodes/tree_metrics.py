"""
Williamson & Shmoys, Chapter 8.5: Probabilistic Approximation of Metrics by Tree Metrics
Fakcharoenphol-Rao-Talwar (FRT) theorem: embed any n-point metric into a
distribution over tree metrics with expected distortion O(log n).
"""

import random
import math
from collections import defaultdict


def frt_embedding(points, dist_func, alpha=4.0, max_depth=15):
    """
    FRT embedding using recursive random clustering.

    Returns dict mapping point -> list of (ancestor, weight) forming path to root.
    """
    n = len(points)
    if n == 0:
        return {}
    if n == 1:
        return {points[0]: []}

    def _build(pts, depth=0):
        if len(pts) <= 1 or depth >= max_depth:
            return {p: [] for p in pts}

        max_d = max(dist_func(pts[i], pts[j])
                     for i in range(len(pts)) for j in range(i+1, len(pts)))
        if max_d == 0:
            return {p: [] for p in pts}

        r = max_d / alpha

        centers = [pts[0]]
        for p in pts[1:]:
            if all(dist_func(p, c) >= r for c in centers):
                centers.append(p)

        clusters = defaultdict(list)
        for p in pts:
            best = min(centers, key=lambda c: dist_func(p, c))
            clusters[best].append(p)

        subtree_paths = {}
        for c in centers:
            cp = clusters[c]
            sub = _build(cp, depth + 1)
            for p in cp:
                subtree_paths[p] = sub[p]

        result = {}
        for c in centers:
            for p in clusters[c]:
                d = dist_func(p, c)
                rd = r * math.ceil(d / r) if r > 0 else d
                path_to_c = subtree_paths.get(p, [])
                result[p] = [(c, rd)] + path_to_c
                if p == c:
                    result[p] = subtree_paths.get(p, [])

        return result

    paths = _build(points)

    all_center_ids = set()
    def collect_centers(pts, depth=0):
        if len(pts) <= 1 or depth >= max_depth:
            return
        max_d = max(dist_func(pts[i], pts[j])
                     for i in range(len(pts)) for j in range(i+1, len(pts)))
        if max_d == 0:
            return
        r = max_d / alpha
        centers = [pts[0]]
        for p in pts[1:]:
            if all(dist_func(p, c) >= r for c in centers):
                centers.append(p)
        for c in centers:
            all_center_ids.add(c)
        clusters = defaultdict(list)
        for p in pts:
            best = min(centers, key=lambda c: dist_func(p, c))
            clusters[best].append(p)
        for c in centers:
            collect_centers(clusters[c], depth + 1)

    collect_centers(points)

    tree = {}
    for p, path in paths.items():
        if path:
            tree[p] = (path[0][0], path[0][1])
        else:
            tree[p] = (None, 0)

    top_centers = [p for p in all_center_ids if tree.get(p, (None,))[0] is None]
    if len(top_centers) > 1:
        virtual_root = -999
        tree[virtual_root] = (None, 0)
        for c in top_centers:
            tree[c] = (virtual_root, 0)

    return tree


def compute_tree_distance(tree, u, v):
    """Compute distance between u and v in the tree metric."""
    if u == v:
        return 0

    def get_path(x):
        path = [x]
        while x in tree and tree[x][0] is not None:
            x = tree[x][0]
            path.append(x)
        return path

    pu = get_path(u)
    pv = get_path(v)
    set_v = set(pv)

    lca = None
    for p in pu:
        if p in set_v:
            lca = p
            break

    if lca is None:
        return float('inf')

    d1 = 0
    x = u
    while x != lca and x in tree and tree[x][0] is not None:
        d1 += tree[x][1]
        x = tree[x][0]

    d2 = 0
    x = v
    while x != lca and x in tree and tree[x][0] is not None:
        d2 += tree[x][1]
        x = tree[x][0]

    return d1 + d2


def evaluate_embedding(points, dist_func, num_samples=10):
    max_ratio = 0
    total_ratio = 0
    count = 0

    for _ in range(num_samples):
        tree = frt_embedding(points, dist_func)
        for i in range(len(points)):
            for j in range(i+1, len(points)):
                orig = dist_func(points[i], points[j])
                if orig == 0:
                    continue
                td = compute_tree_distance(tree, points[i], points[j])
                if td < float('inf') and td > 0:
                    ratio = td / orig
                    total_ratio += ratio
                    count += 1
                    if ratio > max_ratio:
                        max_ratio = ratio

    return max_ratio, total_ratio / count if count > 0 else 0


def demo():
    print("=" * 60)
    print("FRT Tree Metric Embedding (Williamson & Shmoys Ch 8.5)")
    print("Expected distortion O(log n)")
    print("=" * 60)

    points = [0, 1, 2, 3, 4, 5]
    dist_matrix = {
        (0, 1): 1, (0, 2): 3, (0, 3): 2, (0, 4): 4, (0, 5): 5,
        (1, 2): 2, (1, 3): 3, (1, 4): 5, (1, 5): 6,
        (2, 3): 1, (2, 4): 3, (2, 5): 4,
        (3, 4): 2, (3, 5): 3,
        (4, 5): 1,
    }

    def dist_func(u, v):
        if u == v: return 0
        return dist_matrix[(min(u, v), max(u, v))]

    print(f"\nMetric space: {len(points)} points")
    print(f"Distance matrix:")
    for i in range(len(points)):
        print(f"  {[dist_func(points[i], points[j]) for j in range(len(points))]}")

    random.seed(42)
    tree = frt_embedding(points, dist_func)
    print(f"\nSample tree:")
    for p in sorted(tree.keys()):
        parent, weight = tree[p]
        print(f"  {p} -> {parent} (weight {weight:.1f})")

    print(f"\nTree distances:")
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            orig = dist_func(points[i], points[j])
            td = compute_tree_distance(tree, points[i], points[j])
            r = f"{td/orig:.2f}" if orig > 0 and td < float('inf') else "inf"
            print(f"  d({points[i]},{points[j]}): orig={orig}, tree={td:.1f}, ratio={r}")

    max_r, avg_r = evaluate_embedding(points, dist_func, num_samples=10)
    print(f"\n--- Distortion (10 samples) ---")
    print(f"Max ratio: {max_r:.2f}")
    print(f"Average ratio: {avg_r:.2f}")
    print(f"O(log n) bound: {math.log(len(points)):.2f}")


if __name__ == "__main__":
    demo()

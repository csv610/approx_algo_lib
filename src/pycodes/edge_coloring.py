"""
Williamson & Shmoys, Chapter 2.7: Edge Coloring
Exact algorithm: finds a (∆+1)-edge-coloring using fan sequences and local recoloring.

Given graph G = (V, E), color edges so no two edges sharing an endpoint get the same color.
Uses at most ∆+1 colors where ∆ = max degree.
"""

from collections import defaultdict


def edge_coloring(n, edges):
    """
    Edge coloring using the greedy + fan-sequence recoloring algorithm.

    Args:
        n: number of vertices (1-indexed: vertices 1..n)
        edges: list of (u, v) pairs

    Returns:
        dict mapping (u, v) -> color (1-indexed), using at most Δ+1 colors
    """
    adj = defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)

    delta = max(len(adj[v]) for v in range(1, n + 1)) if n > 0 else 0
    num_colors = delta + 1

    color_of = {}
    edge_color = {}
    for v in range(1, n + 1):
        color_of[v] = set()

    def get_available_color(v):
        for c in range(1, num_colors + 1):
            if c not in color_of[v]:
                return c
        return None

    def colors_lacking(v):
        return [c for c in range(1, num_colors + 1) if c not in color_of[v]]

    def get_neighbor_with_color(u, c):
        for w in adj[u]:
            if edge_color.get((min(u, w), max(u, w))) == c:
                return w
        return None

    def switch_colors(v, c1, c2):
        visited = {v}
        stack = [v]
        component = []
        while stack:
            x = stack.pop()
            component.append(x)
            for w in adj[x]:
                if w not in visited:
                    ec = edge_color.get((min(x, w), max(x, w)))
                    if ec == c1 or ec == c2:
                        visited.add(w)
                        stack.append(w)
        for x in component:
            for w in adj[x]:
                me = (min(x, w), max(x, w))
                if me in edge_color:
                    ec = edge_color[me]
                    if ec == c1:
                        edge_color[me] = c2
                        color_of[x].discard(c1)
                        color_of[x].add(c2)
                        color_of[w].discard(c1)
                        color_of[w].add(c2)
                    elif ec == c2:
                        edge_color[me] = c1
                        color_of[x].discard(c2)
                        color_of[x].add(c1)
                        color_of[w].discard(c2)
                        color_of[w].add(c1)

    def color_edge(u, v, c):
        me = (min(u, v), max(u, v))
        edge_color[me] = c
        color_of[u].add(c)
        color_of[v].add(c)

    def build_fan(u, v0):
        """Build fan sequence starting from (u, v0). Returns (fan, colors_used)."""
        fan = [v0]
        colors_used = []

        while True:
            vi = fan[-1]
            lacking_u = colors_lacking(u)
            lacking_vi = colors_lacking(vi)

            common = [c for c in lacking_u if c in lacking_vi]
            if common:
                ci = common[0]
                colors_used.append(ci)
                return fan, colors_used, ci, True

            ci = lacking_vi[0] if lacking_vi else None
            if ci is None:
                break
            colors_used.append(ci)
            wi = get_neighbor_with_color(u, ci)
            if wi is None:
                break
            if wi in fan:
                idx = fan.index(wi)
                return fan, colors_used, ci, False
            fan.append(wi)

        return fan, colors_used, None, None

    uncolored_edges = [(min(u, v), max(u, v)) for u, v in edges]

    for me in uncolored_edges:
        u, v = me
        if me in edge_color:
            continue

        fan, colors_used, ci, case1 = build_fan(u, v)

        if case1 is True:
            color_edge(u, v, ci)
        elif case1 is False:
            wi = fan[-1]
            cu = get_available_color(u)
            if cu is None:
                color_edge(u, v, ci)
                continue

            idx = fan.index(wi)
            if idx > 0:
                for k in range(idx - 1, -1, -1):
                    vk = fan[k]
                    ck = colors_used[k] if k < len(colors_used) else None
                    ck1 = colors_used[k + 1] if k + 1 < len(colors_used) else None
                    if ck1 is not None:
                        me2 = (min(u, vk), max(u, vk))
                        old_color = edge_color.pop(me2, None)
                        if old_color:
                            color_of[u].discard(old_color)
                            color_of[vk].discard(old_color)
                        color_edge(u, vk, ck1)

            switch_colors(u, cu, ci)
            color_edge(u, wi, cu)
        else:
            color_edge(u, v, ci if ci else 1)

    for me in list(edge_color.keys()):
        if edge_color[me] is None:
            u, v = me
            c = get_available_color(u)
            if c is None:
                c = 1
            edge_color[me] = c
            color_of[u].add(c)
            color_of[v].add(c)

    return edge_color


def verify_coloring(n, edges, edge_color):
    """Verify no two adjacent edges share a color."""
    adj = defaultdict(list)
    for u, v in edges:
        me = (min(u, v), max(u, v))
        c = edge_color.get(me)
        for w in (u, v):
            for prev_u, prev_v in adj[w]:
                pm = (min(prev_u, prev_v), max(prev_u, prev_v))
                if edge_color.get(pm) == c:
                    return False
        adj[u].append((u, v))
        adj[v].append((u, v))
    return True


def demo():
    print("=" * 60)
    print("Edge Coloring (Williamson & Shmoys Ch 2.7)")
    print("=" * 60)

    n = 6
    edges = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
    print(f"\nGraph: {n} vertices, {len(edges)} edges")

    adj = defaultdict(int)
    for u, v in edges:
        adj[u] += 1
        adj[v] += 1
    delta = max(adj.values())
    print(f"Max degree Δ = {delta}, using at most {delta + 1} colors")

    coloring = edge_coloring(n, edges)
    print(f"\nEdge coloring:")
    for me in sorted(coloring.keys()):
        print(f"  Edge {me[0]}-{me[1]}: color {coloring[me]}")

    valid = verify_coloring(n, edges, coloring)
    print(f"\nValid coloring: {valid}")
    print(f"Colors used: {len(set(coloring.values()))}")


if __name__ == "__main__":
    demo()

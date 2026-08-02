"""
Chapter 20: Multicut in General Graphs
======================================
Vazirani Ch. 20: Multicut in General Graphs via LP Rounding.

The Multicut Problem:
Given an undirected graph G=(V,E) with edge costs, and k demand pairs
(s_i, t_i), find a minimum cost set of edges whose removal disconnects
all demand pairs. NP-hard in general graphs.

Vazirani's approach:
1. LP relaxation with path constraints
2. Region growing technique for O(log k) approximation

LP Relaxation:
  min  sum c_e * x_e
  s.t. sum_{e in P} x_e >= 1  for all s_i-t_i paths P, all pairs i
       0 <= x_e <= 1

Dual:
  max  sum_j y_j
  s.t. sum_{j: e in path_j} y_j <= c_e  for all edges e
       y_j >= 0

Implements:
1. LP relaxation solver via dual formulation (enumerate paths)
2. Region growing algorithm (O(log k) approximation)
3. Greedy algorithm for comparison
"""

from typing import List, Tuple, Dict, Set
from lp_algorithms import Simplex
import heapq


# ============================================================
# Helpers
# ============================================================

def _build_adj(n: int, edges: List[Tuple[int, int]]) -> Dict[int, List[int]]:
    """Build adjacency list from edge list."""
    adj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def _edge_index(edges: List[Tuple[int, int]]) -> Dict[Tuple[int, int], int]:
    """Map sorted edge tuple to its index."""
    return {tuple(sorted(e)): i for i, e in enumerate(edges)}


def _find_all_simple_paths(
    n: int,
    adj: Dict[int, List[int]],
    src: int,
    dst: int,
    max_paths: int = 200,
) -> List[List[Tuple[int, int]]]:
    """Find all simple paths from src to dst (DFS), returned as edge lists."""
    paths: List[List[Tuple[int, int]]] = []

    def dfs(u: int, visited: Set[int], path: List[Tuple[int, int]]) -> None:
        if len(paths) >= max_paths:
            return
        if u == dst:
            paths.append(list(path))
            return
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                path.append(tuple(sorted((u, v))))
                dfs(v, visited, path)
                path.pop()
                visited.remove(v)

    dfs(src, {src}, [])
    return paths


def _is_connected(
    n: int,
    adj: Dict[int, List[int]],
    s: int,
    t: int,
    removed: Set[Tuple[int, int]],
) -> bool:
    """Check if s and t are connected in G minus the removed edges."""
    visited = {s}
    stack = [s]
    while stack:
        u = stack.pop()
        if u == t:
            return True
        for v in adj[u]:
            e = tuple(sorted((u, v)))
            if e not in removed and v not in visited:
                visited.add(v)
                stack.append(v)
    return False


# ============================================================
# LP Relaxation
# ============================================================

def solve_multicut_lp(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    pairs: List[Tuple[int, int]],
) -> Tuple[List[float], float]:
    """
    Solve LP relaxation for multicut via dual formulation.

    Primal:
      min  sum c_e * x_e
      s.t. sum_{e in P} x_e >= 1   (for every s_i-t_i path P)
           x_e >= 0

    Dual:
      max  sum_j y_j
      s.t. sum_{j: e in path_j} y_j <= c_e   (for every edge e)
           y_j >= 0

    We enumerate simple s_i-t_i paths (bounded by max_paths) and solve
    the dual with Simplex. Primal x_e values are extracted from the
    reduced costs of the slack variables in the final tableau.

    Args:
        n: Number of vertices (0..n-1).
        edges: List of undirected edges.
        costs: Edge costs (same order as edges).
        pairs: Demand pairs (s_i, t_i).

    Returns:
        (x_values, lp_optimal_value) where x_e in [0, 1].
    """
    m = len(edges)
    k = len(pairs)

    if k == 0:
        return [0.0] * m, 0.0

    adj = _build_adj(n, edges)
    eidx = _edge_index(edges)

    # Enumerate all simple paths for each demand pair
    path_list: List[Tuple[int, List[Tuple[int, int]]]] = []  # (pair_idx, path_edges)
    for i, (s, t) in enumerate(pairs):
        paths = _find_all_simple_paths(n, adj, s, t)
        for p in paths:
            path_list.append((i, p))

    if not path_list:
        # No paths exist; pairs are already disconnected
        return [0.0] * m, 0.0

    total_paths = len(path_list)

    # Build dual constraint matrix  A_dual : m x total_paths
    # A_dual[e][j] = 1  iff  edge e appears in path j
    A_dual: List[List[float]] = [[0.0] * total_paths for _ in range(m)]
    for col_j, (_pair_idx, path) in enumerate(path_list):
        for edge in path:
            if edge in eidx:
                A_dual[eidx[edge]][col_j] = 1.0

    b_dual = costs[:]                     # c_e for each edge
    c_dual = [1.0] * total_paths          # maximize sum y_j

    simplex = Simplex(A_dual, b_dual, c_dual)
    y, dual_opt = simplex.solve()

    if y is None:
        # LP infeasible or unbounded (should not happen for valid instances)
        return [0.0] * m, float("inf")

    # Extract primal x from slack-variable reduced costs in the final tableau
    x = simplex.obj_row[total_paths : total_paths + m]

    # Clamp to valid range [0, 1]
    x = [max(0.0, min(1.0, xi)) for xi in x]

    return x, dual_opt


# ============================================================
# Region Growing  (O(log k) approximation)
# ============================================================

def multicut_region_growing(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    pairs: List[Tuple[int, int]],
) -> Tuple[List[Tuple[int, int]], float]:
    """
    Multicut via region growing (Vazirani Ch. 20).

    Algorithm:
    1. Solve the LP relaxation to obtain fractional edge values x*.
    2. For each demand pair (s_i, t_i):
       a. Grow a ball B_i from s_i using Dijkstra with edge weights
          w_e = c_e / x*_e  (edges with x*_e ~ 0 are effectively infinite).
       b. Stop the ball the moment t_i is reached.
       c. Every edge crossing the boundary (B_i, V \\ B_i) is added to the cut.
    3. Return the union of all boundary edges.

    The scaling argument: by LP feasibility every s_i-t_i path has
    sum x*_e >= 1, so the ball must cross at least one edge with
    sufficient x* weight.  Repeating with a random ordering of the
    pairs yields an expected cost of O(ln k) * LP_OPT.

    Args:
        n: Number of vertices (0..n-1).
        edges: List of undirected edges.
        costs: Edge costs (same order as edges).
        pairs: Demand pairs (s_i, t_i).

    Returns:
        (cut_edges, total_cost)
    """
    m = len(edges)
    k = len(pairs)

    if k == 0:
        return [], 0.0

    # Step 1: Solve LP
    x, lp_val = solve_multicut_lp(n, edges, costs, pairs)

    adj = _build_adj(n, edges)
    eidx = _edge_index(edges)

    # Step 2: Region growing for each pair
    chosen_edges: Set[Tuple[int, int]] = set()
    INF = float("inf")

    for s, t in pairs:
        # Dijkstra from s with weights c_e / x_e
        dist = [INF] * n
        dist[s] = 0.0
        visited = [False] * n
        pq: List[Tuple[float, int]] = [(0.0, s)]

        while pq:
            d, u = heapq.heappop(pq)
            if visited[u]:
                continue
            visited[u] = True
            if u == t:
                break
            for v in adj[u]:
                if not visited[v]:
                    e = tuple(sorted((u, v)))
                    e_idx = eidx[e]
                    xe = max(x[e_idx], 1e-9)
                    w = costs[e_idx] / xe
                    nd = d + w
                    if nd < dist[v]:
                        dist[v] = nd
                        heapq.heappush(pq, (nd, v))

        # Ball B_i = {v : dist[v] < dist[t]}  (strictly before t)
        ball = {v for v in range(n) if dist[v] < dist[t] - 1e-9}

        # Add boundary edges
        for u_val, v_val in edges:
            if (u_val in ball) != (v_val in ball):
                chosen_edges.add(tuple(sorted((u_val, v_val))))

    cut_edges = list(chosen_edges)
    total_cost = sum(costs[eidx[e]] for e in cut_edges)
    return cut_edges, total_cost


# ============================================================
# Greedy  (comparison heuristic)
# ============================================================

def multicut_greedy(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    pairs: List[Tuple[int, int]],
) -> Tuple[List[Tuple[int, int]], float]:
    """
    Greedy heuristic for multicut.

    Algorithm:
    1. Sort edges by cost (ascending).
    2. Add edges one by one; after each addition, check whether all
       demand pairs are disconnected.
    3. Stop as soon as every pair is cut.

    In highly connected graphs no single edge removal may disconnect
    any pair, so the greedy simply accumulates cheap edges until the
    cut property is satisfied.  It has no worst-case guarantee but
    provides a baseline for comparison.

    Args:
        n: Number of vertices (0..n-1).
        edges: List of undirected edges.
        costs: Edge costs (same order as edges).
        pairs: Demand pairs (s_i, t_i).

    Returns:
        (cut_edges, total_cost)
    """
    if not pairs:
        return [], 0.0

    adj = _build_adj(n, edges)
    eidx = _edge_index(edges)

    # Sort edge indices by cost (ascending)
    sorted_indices = sorted(range(len(edges)), key=lambda i: costs[i])

    cut_set: Set[Tuple[int, int]] = set()

    for idx in sorted_indices:
        e = tuple(sorted(edges[idx]))
        cut_set.add(e)

        # Check if all pairs are now disconnected
        if all(
            not _is_connected(n, adj, s, t, cut_set) for s, t in pairs
        ):
            break

    cut_edges = list(cut_set)
    total_cost = sum(costs[eidx[e]] for e in cut_edges)
    return cut_edges, total_cost


# ============================================================
# Verification helper
# ============================================================

def _verify_cut(
    n: int,
    edges: List[Tuple[int, int]],
    cut_edges: List[Tuple[int, int]],
    pairs: List[Tuple[int, int]],
) -> bool:
    """Return True if every demand pair is disconnected by the cut."""
    adj = _build_adj(n, edges)
    cut_set = set(tuple(sorted(e)) for e in cut_edges)
    for s, t in pairs:
        if _is_connected(n, adj, s, t, cut_set):
            return False
    return True


# ============================================================
# Demo
# ============================================================

def demo_multicut_general() -> None:
    """Demo for Multicut in General Graphs (Vazirani Ch. 20)."""
    print("=" * 60)
    print("Chapter 20: Multicut in General Graphs")
    print("=" * 60)

    # ---- Instance 1: diamond graph with cross edges ----
    #    0 ---1--- 1
    #    |         |
    #    2    X    2
    #    |         |
    #    2 ---1--- 3
    n = 4
    edges = [
        (0, 1),  # cost 1
        (0, 2),  # cost 2
        (1, 3),  # cost 2
        (2, 3),  # cost 1
        (0, 3),  # cost 5  (cross)
        (1, 2),  # cost 5  (cross)
    ]
    costs = [1.0, 2.0, 2.0, 1.0, 5.0, 5.0]
    pairs = [(0, 3), (1, 2)]

    print("\n--- Instance 1: Diamond graph ---")
    print(f"  Vertices: {n}")
    print("  Edges & Costs:")
    for i, (u, v) in enumerate(edges):
        print(f"    ({u}, {v}) : cost = {costs[i]:.1f}")
    print(f"  Demand Pairs: {pairs}")

    # LP
    x, lp_val = solve_multicut_lp(n, edges, costs, pairs)
    print(f"\n1. LP Relaxation:")
    print(f"  LP Optimal Value: {lp_val:.4f}")
    print(f"  Edge values x*:")
    for i, (u, v) in enumerate(edges):
        print(f"    ({u},{v}) : x* = {x[i]:.4f}")

    # Region Growing
    cut_rg, cost_rg = multicut_region_growing(n, edges, costs, pairs)
    valid_rg = _verify_cut(n, edges, cut_rg, pairs)
    print(f"\n2. Region Growing:")
    print(f"  Cut edges: {cut_rg}")
    print(f"  Total cost: {cost_rg:.2f}")
    print(f"  Valid cut: {valid_rg}")

    # Greedy
    cut_gr, cost_gr = multicut_greedy(n, edges, costs, pairs)
    valid_gr = _verify_cut(n, edges, cut_gr, pairs)
    print(f"\n3. Greedy:")
    print(f"  Cut edges: {cut_gr}")
    print(f"  Total cost: {cost_gr:.2f}")
    print(f"  Valid cut: {valid_gr}")

    # Comparison
    print(f"\n4. Comparison:")
    print(f"  LP Lower Bound: {lp_val:.4f}")
    if lp_val > 0:
        print(f"  Region Growing: {cost_rg:.2f}  (approx ratio: {cost_rg / lp_val:.2f}x)")
        print(f"  Greedy:         {cost_gr:.2f}  (approx ratio: {cost_gr / lp_val:.2f}x)")
    else:
        print("  (LP bound is zero; pairs may already be disconnected)")

    # ---- Instance 2: cycle of 6 ----
    print("\n" + "=" * 60)
    print("--- Instance 2: Cycle C6 ---")
    n2 = 6
    edges2 = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
    ]
    costs2 = [1.0, 3.0, 1.0, 3.0, 1.0, 3.0]
    pairs2 = [(0, 3), (1, 4), (2, 5)]

    print(f"  Vertices: {n2}")
    print("  Edges & Costs:")
    for i, (u, v) in enumerate(edges2):
        print(f"    ({u}, {v}) : cost = {costs2[i]:.1f}")
    print(f"  Demand Pairs: {pairs2}")

    x2, lp2 = solve_multicut_lp(n2, edges2, costs2, pairs2)
    print(f"\n  LP Optimal: {lp2:.4f}")
    print(f"  Edge values: {['%.3f' % v for v in x2]}")

    cut_rg2, cost_rg2 = multicut_region_growing(n2, edges2, costs2, pairs2)
    valid_rg2 = _verify_cut(n2, edges2, cut_rg2, pairs2)
    print(f"\n  Region Growing: cost = {cost_rg2:.2f}, valid = {valid_rg2}")
    print(f"    Cut: {cut_rg2}")

    cut_gr2, cost_gr2 = multicut_greedy(n2, edges2, costs2, pairs2)
    valid_gr2 = _verify_cut(n2, edges2, cut_gr2, pairs2)
    print(f"  Greedy:         cost = {cost_gr2:.2f}, valid = {valid_gr2}")
    print(f"    Cut: {cut_gr2}")

    if lp2 > 0:
        print(f"\n  Ratios: RG={cost_rg2 / lp2:.2f}x  Greedy={cost_gr2 / lp2:.2f}x")

    # ---- Instance 3: K4 ----
    print("\n" + "=" * 60)
    print("--- Instance 3: Complete graph K4 ---")
    n3 = 4
    edges3 = [
        (0, 1), (0, 2), (0, 3),
        (1, 2), (1, 3), (2, 3),
    ]
    costs3 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    pairs3 = [(0, 1), (2, 3)]

    print(f"  Vertices: {n3}")
    print("  Edges & Costs:")
    for i, (u, v) in enumerate(edges3):
        print(f"    ({u}, {v}) : cost = {costs3[i]:.1f}")
    print(f"  Demand Pairs: {pairs3}")

    x3, lp3 = solve_multicut_lp(n3, edges3, costs3, pairs3)
    print(f"\n  LP Optimal: {lp3:.4f}")

    cut_rg3, cost_rg3 = multicut_region_growing(n3, edges3, costs3, pairs3)
    valid_rg3 = _verify_cut(n3, edges3, cut_rg3, pairs3)
    print(f"  Region Growing: cost = {cost_rg3:.2f}, valid = {valid_rg3}")

    cut_gr3, cost_gr3 = multicut_greedy(n3, edges3, costs3, pairs3)
    valid_gr3 = _verify_cut(n3, edges3, cut_gr3, pairs3)
    print(f"  Greedy:         cost = {cost_gr3:.2f}, valid = {valid_gr3}")

    if lp3 > 0:
        print(f"\n  Ratios: RG={cost_rg3 / lp3:.2f}x  Greedy={cost_gr3 / lp3:.2f}x")

    print()


if __name__ == "__main__":
    demo_multicut_general()

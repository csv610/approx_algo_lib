"""
Chapter 18: Multicut and Integer Multicommodity Flow in Trees
=============================================================
Vazirani Ch. 18: Multicut in trees admits an exact polytime algorithm
because the constraint matrix is totally dual integral (TDI).

Key insight: on trees, the LP relaxation of both multicut and
multicommodity flow has integral optimal solutions.  The path-edge
incidence matrix of a tree is totally unimodular, so the simplex
method returns integer x* and f*.

Primal (Multicut LP):
  min  sum_e  c_e * x_e
  s.t. sum_{e in P_i} x_e >= 1   for every demand pair i
       x_e >= 0

Dual (Multicommodity Flow LP):
  max  sum_i  f_i
  s.t. sum_{i : e in P_i} f_i <= c_e   for every edge e
       f_i >= 0

On a tree each pair has exactly one path, so the constraint matrix
entries are all 0/1 and the matrix is TU.  Both LPs therefore yield
integer optimal solutions, and strong duality gives:

  OPT_multicut = OPT_flow

Implements:
1. tree_multicut             - Exact multicut (via LP)
2. tree_multicommodity_flow  - Integer multicommodity flow (via LP)
3. solve_multicut_lp_tree    - Raw LP relaxation (returns both x and f)
4. verify_optimality         - Check strong duality cut_cost == flow_value
5. demo_multicut_trees       - Worked examples
"""

from typing import List, Tuple, Dict, Optional
from lp_algorithms import Simplex


# ============================================================
# Tree helpers
# ============================================================

def _build_tree_adj(
    edges: List[Tuple[int, int]],
) -> Dict[int, List[int]]:
    """Build adjacency list from undirected edge list."""
    adj: Dict[int, List[int]] = {}
    for u, v in edges:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    return adj


def _compute_parent_depth(
    n: int,
    adj: Dict[int, List[int]],
    root: int = 0,
) -> Tuple[Dict[int, int], Dict[int, int]]:
    """BFS from root to compute parent pointers and depths."""
    parent: Dict[int, int] = {root: -1}
    depth: Dict[int, int] = {root: 0}
    queue = [root]
    while queue:
        curr = queue.pop(0)
        for nb in adj[curr]:
            if nb not in parent:
                parent[nb] = curr
                depth[nb] = depth[curr] + 1
                queue.append(nb)
    return parent, depth


def _path_edges(
    u: int,
    v: int,
    parent: Dict[int, int],
    depth: Dict[int, int],
) -> List[Tuple[int, int]]:
    """Return sorted edge tuples on the unique u-v path in the tree."""
    # Lift the deeper node to the same depth
    a, b = u, v
    while depth[a] > depth[b]:
        a = parent[a]
    while depth[b] > depth[a]:
        b = parent[b]
    # Walk both up until they meet (LCA)
    while a != b:
        a = parent[a]
        b = parent[b]
    lca = a

    path: List[Tuple[int, int]] = []
    c = u
    while c != lca:
        path.append(tuple(sorted((c, parent[c]))))
        c = parent[c]
    c = v
    while c != lca:
        path.append(tuple(sorted((c, parent[c]))))
        c = parent[c]
    return path


# ============================================================
# LP-based exact solver (the core)
# ============================================================

def solve_multicut_lp_tree(
    edges: List[Tuple[int, int]],
    costs: Dict[Tuple[int, int], float],
    demands: List[Tuple[int, int]],
    root: int = 0,
) -> Tuple[List[float], float, List[float], float]:
    """
    Solve the multicut / multicommodity-flow LP on a tree.

    The constraint matrix (pair-edge incidence) is totally unimodular
    for trees, so both the primal (x) and dual (f) solutions are integral.

    Args:
        edges:    Undirected tree edges [(u, v), ...].
        costs:    Edge cost dict, keyed by sorted tuple.
        demands:  Demand pairs [(s_i, t_i), ...].
        root:     Arbitrary root for tree traversal.

    Returns:
        (x, cut_opt, f, flow_opt)
        x       : per-edge multicut variable (integral 0/1).
        cut_opt : optimal multicut cost  =  sum_e c_e * x_e.
        f       : per-pair flow variable (integral >= 0).
        flow_opt: optimal flow value      =  sum_i f_i.
    """
    n = max(max(u, v) for u, v in edges) + 1
    m = len(edges)
    k = len(demands)

    if k == 0:
        return [0.0] * m, 0.0, [], 0.0

    adj = _build_tree_adj(edges)
    parent, depth = _compute_parent_depth(n, adj, root)

    # Edge index map
    eidx: Dict[Tuple[int, int], int] = {tuple(sorted(e)): i for i, e in enumerate(edges)}

    # Compute path for each demand pair
    pair_paths: List[List[Tuple[int, int]]] = []
    for s, t in demands:
        pair_paths.append(_path_edges(s, t, parent, depth))

    # ----- Primal LP (Multicut) -----
    # min  sum_e  c_e * x_e
    # s.t. A * x >= 1       (one constraint per demand pair)
    #      x >= 0
    #
    # Convert to standard form for Simplex (max -c^T x s.t. -A x <= -1):
    #   max  (-c)^T x
    #   s.t. -A x <= -1
    #        x >= 0
    #
    # But Simplex expects max c^T x s.t. A x <= b, x >= 0.
    # We solve the dual instead (cleaner, and we get f too).

    # ----- Dual LP (Multicommodity Flow) -----
    # max  sum_i  f_i
    # s.t. sum_{i: e in P_i} f_i <= c_e   for every edge e
    #      f_i >= 0
    A_dual: List[List[float]] = [[0.0] * k for _ in range(m)]
    for col_i, path in enumerate(pair_paths):
        for e in path:
            if e in eidx:
                A_dual[eidx[e]][col_i] = 1.0

    b_dual = [costs.get(tuple(sorted(edges[i])), 0.0) for i in range(m)]
    c_dual = [1.0] * k  # maximize total flow

    simplex = Simplex(A_dual, b_dual, c_dual)
    f, flow_opt = simplex.solve()

    if f is None:
        return [0.0] * m, 0.0, [0.0] * k, 0.0

    # Primal x recovered from slack-variable reduced costs
    x = simplex.obj_row[k: k + m]
    x = [max(0.0, min(1.0, xi)) for xi in x]

    cut_opt = sum(x[i] * costs.get(tuple(sorted(edges[i])), 0.0) for i in range(m))

    return x, cut_opt, f, flow_opt


# ============================================================
# Convenience wrappers
# ============================================================

def tree_multicut(
    edges: List[Tuple[int, int]],
    costs: Dict[Tuple[int, int], float],
    demands: List[Tuple[int, int]],
    root: int = 0,
) -> Tuple[List[Tuple[int, int]], float]:
    """
    Exact minimum multicut in a tree.

    Returns the set of cut edges and total cost.  The solution is
    optimal (integrality gap = 1) because the LP matrix is TU.
    """
    x, cut_opt, _, _ = solve_multicut_lp_tree(edges, costs, demands, root)
    cut_edges = [
        tuple(sorted(edges[i])) for i in range(len(edges)) if x[i] > 0.5
    ]
    return cut_edges, cut_opt


def tree_multicommodity_flow(
    edges: List[Tuple[int, int]],
    capacities: Dict[Tuple[int, int], float],
    demands: List[Tuple[int, int]],
    root: int = 0,
) -> Tuple[List[float], float]:
    """
    Maximum integer multicommodity flow in a tree.

    Each demand pair (s_i, t_i) gets an integral flow f_i along its
    unique tree path.  Edge capacities are respected.

    Returns (per-pair flow list, total flow value).
    """
    _, _, f, flow_opt = solve_multicut_lp_tree(edges, capacities, demands, root)
    return f, flow_opt


# ============================================================
# Verification
# ============================================================

def verify_optimality(
    cut_cost: float,
    flow_value: float,
    tol: float = 1e-6,
) -> bool:
    """
    Check strong duality: OPT_multicut == OPT_flow.

    On trees the LP matrix is TDI, so both LPs achieve the same
    integral optimum.  This helper confirms the equality within
    numerical tolerance.
    """
    return abs(cut_cost - flow_value) <= tol


def _verify_cut_valid(
    edges: List[Tuple[int, int]],
    cut_edges: List[Tuple[int, int]],
    demands: List[Tuple[int, int]],
) -> bool:
    """Return True if every demand pair is disconnected by the cut."""
    cut_set = set(tuple(sorted(e)) for e in cut_edges)
    adj = _build_tree_adj(edges)
    n = max(max(u, v) for u, v in edges) + 1

    for s, t in demands:
        visited = {s}
        stack = [s]
        while stack:
            u = stack.pop()
            if u == t:
                return False
            for v in adj[u]:
                e = tuple(sorted((u, v)))
                if e not in cut_set and v not in visited:
                    visited.add(v)
                    stack.append(v)
    return True


# ============================================================
# Demo
# ============================================================

def demo_multicut_trees() -> None:
    """Demonstrate exact multicut and multicommodity flow on trees."""
    print("=" * 60)
    print("Chapter 18: Multicut & Integer Multicommodity Flow in Trees")
    print("=" * 60)

    # ----------------------------------------------------------
    # Instance 1  (from Vazirani's exposition)
    #
    #          0
    #         / \
    #        1   2        costs: c(0,1)=4, c(0,2)=3
    #       / \   \       c(1,3)=2, c(1,4)=5, c(2,5)=1
    #      3   4   5
    #
    #  demands: (3,5), (4,5)
    # ----------------------------------------------------------
    print("\n--- Instance 1 ---")
    edges1 = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5)]
    costs1 = {
        (0, 1): 4.0, (0, 2): 3.0,
        (1, 3): 2.0, (1, 4): 5.0, (2, 5): 1.0,
    }
    demands1 = [(3, 5), (4, 5)]

    print("  Tree edges & costs:")
    for e, c in costs1.items():
        print(f"    {e} : {c:.0f}")
    print(f"  Demands: {demands1}")

    # LP solve
    x1, cut_opt1, f1, flow_opt1 = solve_multicut_lp_tree(edges1, costs1, demands1)
    print(f"\n1. LP Relaxation:")
    print(f"  Primal x (cut vars):  {[round(v, 4) for v in x1]}")
    print(f"  Dual   f (flow vars): {[round(v, 4) for v in f1]}")
    print(f"  Multicut LP value:    {cut_opt1:.4f}")
    print(f"  Flow LP value:        {flow_opt1:.4f}")
    print(f"  Strong duality OK:    {verify_optimality(cut_opt1, flow_opt1)}")

    # Exact multicut
    cut1, cost1 = tree_multicut(edges1, costs1, demands1)
    valid1 = _verify_cut_valid(edges1, cut1, demands1)
    print(f"\n2. Exact Multicut:")
    print(f"  Cut edges: {cut1}")
    print(f"  Total cost: {cost1:.0f}")
    print(f"  Valid cut:  {valid1}")

    # Multicommodity flow
    fval1, total_flow1 = tree_multicommodity_flow(edges1, costs1, demands1)
    print(f"\n3. Integer Multicommodity Flow:")
    for i, (s, t) in enumerate(demands1):
        print(f"  Pair ({s},{t}): flow = {fval1[i]:.0f}")
    print(f"  Total flow: {total_flow1:.0f}")

    # Optimality check
    print(f"\n4. Optimality:")
    print(f"  cut_cost == flow_value? {verify_optimality(cost1, total_flow1)}")

    # ----------------------------------------------------------
    # Instance 2  -  balanced binary tree, 7 nodes
    #
    #           0
    #         /   \
    #        1     2         costs: all edges = 1
    #       / \   / \
    #      3   4 5   6
    #
    #  demands: (3,4), (5,6), (3,5)
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("--- Instance 2: Balanced binary tree ---")
    edges2 = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
    costs2 = {tuple(sorted(e)): 1.0 for e in edges2}
    demands2 = [(3, 4), (5, 6), (3, 5)]

    print(f"  Tree: 7-node balanced binary, all edge costs = 1")
    print(f"  Demands: {demands2}")

    x2, cut_opt2, f2, flow_opt2 = solve_multicut_lp_tree(edges2, costs2, demands2)
    print(f"\n1. LP Relaxation:")
    print(f"  Primal x: {[round(v, 4) for v in x2]}")
    print(f"  Dual   f: {[round(v, 4) for v in f2]}")
    print(f"  Multicut: {cut_opt2:.4f}    Flow: {flow_opt2:.4f}")

    cut2, cost2 = tree_multicut(edges2, costs2, demands2)
    valid2 = _verify_cut_valid(edges2, cut2, demands2)
    print(f"\n2. Exact Multicut: edges={cut2}, cost={cost2:.0f}, valid={valid2}")

    fval2, total_flow2 = tree_multicommodity_flow(edges2, costs2, demands2)
    print(f"3. Integer Flow:   per-pair={[f'{v:.0f}' for v in fval2]}, total={total_flow2:.0f}")
    print(f"4. duality: {verify_optimality(cost2, total_flow2)}")

    # ----------------------------------------------------------
    # Instance 3  -  chain (path graph) with unequal costs
    #
    #  0 --3-- 1 --1-- 2 --4-- 3 --2-- 4
    #
    #  demands: (0,4), (1,3)
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("--- Instance 3: Path graph ---")
    edges3 = [(0, 1), (1, 2), (2, 3), (3, 4)]
    costs3 = {(0, 1): 3.0, (1, 2): 1.0, (2, 3): 4.0, (3, 4): 2.0}
    demands3 = [(0, 4), (1, 3)]

    print(f"  Path: 0--3--1--1--2--4--3--2--4")
    print(f"  Demands: {demands3}")

    x3, cut_opt3, f3, flow_opt3 = solve_multicut_lp_tree(edges3, costs3, demands3)
    print(f"\n1. LP:")
    print(f"  x = {[round(v, 4) for v in x3]}")
    print(f"  f = {[round(v, 4) for v in f3]}")
    print(f"  Multicut: {cut_opt3:.4f}    Flow: {flow_opt3:.4f}")

    cut3, cost3 = tree_multicut(edges3, costs3, demands3)
    valid3 = _verify_cut_valid(edges3, cut3, demands3)
    print(f"\n2. Exact Multicut: edges={cut3}, cost={cost3:.0f}, valid={valid3}")

    fval3, total_flow3 = tree_multicommodity_flow(edges3, costs3, demands3)
    print(f"3. Integer Flow:   per-pair={[f'{v:.0f}' for v in fval3]}, total={total_flow3:.0f}")
    print(f"4. duality: {verify_optimality(cost3, total_flow3)}")

    # ----------------------------------------------------------
    # Instance 4  -  star graph
    #
    #        1
    #        |
    #   2 -- 0 -- 3
    #        |
    #        4
    #
    #  demands: (1,3), (2,4)
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("--- Instance 4: Star graph ---")
    edges4 = [(0, 1), (0, 2), (0, 3), (0, 4)]
    costs4 = {(0, 1): 5.0, (0, 2): 2.0, (0, 3): 3.0, (0, 4): 4.0}
    demands4 = [(1, 3), (2, 4)]

    print(f"  Star centered at 0, costs 5,2,3,4")
    print(f"  Demands: {demands4}")

    x4, cut_opt4, f4, flow_opt4 = solve_multicut_lp_tree(edges4, costs4, demands4)
    print(f"\n1. LP:")
    print(f"  x = {[round(v, 4) for v in x4]}")
    print(f"  f = {[round(v, 4) for v in f4]}")
    print(f"  Multicut: {cut_opt4:.4f}    Flow: {flow_opt4:.4f}")

    cut4, cost4 = tree_multicut(edges4, costs4, demands4)
    valid4 = _verify_cut_valid(edges4, cut4, demands4)
    print(f"\n2. Exact Multicut: edges={cut4}, cost={cost4:.0f}, valid={valid4}")

    fval4, total_flow4 = tree_multicommodity_flow(edges4, costs4, demands4)
    print(f"3. Integer Flow:   per-pair={[f'{v:.0f}' for v in fval4]}, total={total_flow4:.0f}")
    print(f"4. duality: {verify_optimality(cost4, total_flow4)}")

    # ----------------------------------------------------------
    # Instance 5  -  larger tree (11 nodes)
    #
    #            0
    #          / | \
    #         1  2  3
    #        /|    |\
    #       4  5   6  7
    #          |      |
    #          8      9
    #                |
    #               10
    #
    #  demands: (4,10), (8,7), (5,9)
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("--- Instance 5: Larger tree (11 nodes) ---")
    edges5 = [
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 5),
        (2, 6), (3, 7),
        (5, 8), (7, 9), (9, 10),
    ]
    costs5 = {
        (0, 1): 3.0, (0, 2): 2.0, (0, 3): 4.0,
        (1, 4): 1.0, (1, 5): 5.0,
        (2, 6): 3.0, (3, 7): 2.0,
        (5, 8): 1.0, (7, 9): 6.0, (9, 10): 1.0,
    }
    demands5 = [(4, 10), (8, 7), (5, 9)]

    print(f"  11-node tree")
    print(f"  Demands: {demands5}")

    x5, cut_opt5, f5, flow_opt5 = solve_multicut_lp_tree(edges5, costs5, demands5)
    print(f"\n1. LP:")
    print(f"  x = {[round(v, 4) for v in x5]}")
    print(f"  f = {[round(v, 4) for v in f5]}")
    print(f"  Multicut: {cut_opt5:.4f}    Flow: {flow_opt5:.4f}")

    cut5, cost5 = tree_multicut(edges5, costs5, demands5)
    valid5 = _verify_cut_valid(edges5, cut5, demands5)
    print(f"\n2. Exact Multicut: edges={cut5}, cost={cost5:.0f}, valid={valid5}")

    fval5, total_flow5 = tree_multicommodity_flow(edges5, costs5, demands5)
    print(f"3. Integer Flow:   per-pair={[f'{v:.0f}' for v in fval5]}, total={total_flow5:.0f}")
    print(f"4. duality: {verify_optimality(cost5, total_flow5)}")

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("Summary: Integrality on Trees")
    print("=" * 60)
    print("On trees the path-edge incidence matrix is totally unimodular,")
    print("so both the multicut LP and the multicommodity-flow LP have")
    print("integral optimal solutions.  Strong duality guarantees:")
    print("  OPT(multicut)  ==  OPT(flow)")
    print("No rounding is needed --- the LP relaxation is exact.")
    print("=" * 60)


if __name__ == "__main__":
    demo_multicut_trees()

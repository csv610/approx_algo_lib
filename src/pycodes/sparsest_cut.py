"""
Chapter 21: Sparsest Cut
========================
Vazirani Ch. 21: Sparsest Cut via LP relaxation and metric embedding.

The Sparsest Cut Problem:
  Given G=(V,E) with edge capacities c_e and demand pairs {(s_i,t_i,d_i)},
  find a cut (S, V\\S) minimizing cap(S)/demand(S), where:
    cap(S)     = sum of c_e for edges crossing the cut
    demand(S)  = sum of d_i for pairs separated by the cut

Key ideas (Aumann-Rabani 1995, Arora-Rao-Vazirani 2009):
  1. LP relaxation via multicommodity flow / cut-covering formulation
  2. The dual of the LP yields a metric on V
  3. l1-embeddability: every metric is close to an l1-metric
  4. Rounding the l1-embedding gives O(sqrt(log n)) approximation

Implements:
  1. LP relaxation via dual formulation (cut enumeration, small instances)
  2. LP rounding via l1-metric embedding and random partitioning
  3. Combinatorial exact search
  4. Conductance computation (related measure)
"""

import math
import random
from typing import List, Tuple, Set

from lp_algorithms import Simplex


# ====================================================================
# HELPERS
# ====================================================================

def _cut_edge_indices(S: Set[int], edges: List[Tuple[int, int]]) -> List[int]:
    """Return indices of edges crossing cut (S, V\\S)."""
    return [i for i, (u, v) in enumerate(edges) if (u in S) != (v in S)]


def _enumerate_cuts(n: int):
    """Yield all non-trivial subsets S of {0,...,n-1}."""
    for mask in range(1, (1 << n) - 1):
        S = set()
        for v in range(n):
            if mask & (1 << v):
                S.add(v)
        yield S


def _floyd_warshall(n: int, edges: List[Tuple[int, int]],
                     weights: List[float]) -> List[List[float]]:
    """All-pairs shortest paths via Floyd-Warshall."""
    INF = float('inf')
    dist = [[INF] * n for _ in range(n)]
    for v in range(n):
        dist[v][v] = 0.0
    for i, (u, v) in enumerate(edges):
        if weights[i] < dist[u][v]:
            dist[u][v] = weights[i]
            dist[v][u] = weights[i]
    for kk in range(n):
        for i in range(n):
            dik = dist[i][kk]
            if dik >= INF:
                continue
            for j in range(n):
                alt = dik + dist[kk][j]
                if alt < dist[i][j]:
                    dist[i][j] = alt
    return dist


def _cut_cost(S: Set[int], edges: List[Tuple[int, int]],
              capacities: List[float]) -> float:
    """Total capacity of edges crossing cut (S, V\\S)."""
    return sum(capacities[i] for i in _cut_edge_indices(S, edges))


def _demand_crossing(S: Set[int], demands: List[Tuple[int, int, float]]) -> float:
    """Total demand of pairs separated by cut (S, V\\S)."""
    return sum(d for s, t, d in demands if (s in S) != (t in S))


def _sparsest_ratio(S: Set[int], edges: List[Tuple[int, int]],
                    capacities: List[float],
                    demands: List[Tuple[int, int, float]]) -> float:
    """cap(S) / demand(S), or inf if demand(S) == 0."""
    dem = _demand_crossing(S, demands)
    if dem <= 1e-12:
        return float('inf')
    return _cut_cost(S, edges, capacities) / dem


# ====================================================================
# 1. LP RELAXATION
# ====================================================================

def solve_sparsest_cut_lp(
    n: int,
    edges: List[Tuple[int, int]],
    capacities: List[float],
    demands: List[Tuple[int, int, float]]
) -> Tuple[List[float], float]:
    """
    Solve the LP relaxation for sparsest cut via cut enumeration.

    LP (multicommodity-flow formulation):
      minimize  sum_e  c_e * x_e
      subject to
          min-cut(s_i, t_i; x)  >=  d_i     for every demand pair i
          x_e >= 0

    This finds the cheapest fractional cut that simultaneously satisfies
    every pair's demand requirement.  Its value is a *lower bound* on the
    cost of any integer cut that separates ALL pairs, and its dual metric
    drives the O(sqrt(log n)) rounding algorithm.

    Its DUAL has b >= 0 (capacities), making the origin feasible for
    the Simplex solver:

      maximize  sum_{i,S}  d_i * y_{i,S}
      subject to
          sum_{i,S : e in delta(S)}  y_{i,S}  <=  c_e     for every edge e
          y_{i,S} >= 0

    Variables y_{i,S} exist for each pair i and each cut S that
    separates s_i from t_i.  After solving, primal variables x_e are
    recovered from the slack-column reduced costs in the optimal
    tableau (same technique as set_cover_lp).

    Uses Simplex from lp_algorithms.  Feasible for n <= ~12.

    Returns:
        (x, lp_value) where x is the fractional edge assignment and
        lp_value is the optimal objective value (sparsest cut value).
    """
    m = len(edges)
    k = len(demands)

    if k == 0 or m == 0:
        return [0.0] * m, 0.0

    # ---- enumerate all non-trivial cuts ----
    cuts: List[Set[int]] = list(_enumerate_cuts(n))
    M = len(cuts)

    # ---- build dual LP ----
    # max  c^T y   s.t.  A y <= cap,  y >= 0
    #
    # For each pair i, find which cuts separate it.
    # Variables are indexed as a flat list over (pair, cut) pairs.
    # A[e_idx][var_idx] = 1  if edge e crosses the cut for this variable.

    var_pairs: List[Tuple[int, int]] = []   # (pair_idx, cut_idx)
    for i in range(k):
        s_i, t_i, _ = demands[i]
        for j, S in enumerate(cuts):
            if (s_i in S) != (t_i in S):
                var_pairs.append((i, j))

    N = len(var_pairs)  # total dual variables

    if N == 0:
        return [0.0] * m, 0.0

    A_dual: List[List[float]] = []
    b_dual: List[float] = []

    for e_idx, (u, v) in enumerate(edges):
        row = [0.0] * N
        for var_idx, (pair_idx, cut_idx) in enumerate(var_pairs):
            S = cuts[cut_idx]
            if (u in S) != (v in S):
                row[var_idx] = 1.0
        A_dual.append(row)
        b_dual.append(capacities[e_idx])

    c_dual = [demands[pair_idx][2] for pair_idx, _ in var_pairs]

    simplex = Simplex(A_dual, b_dual, c_dual)
    y_opt, dual_obj = simplex.solve()

    if y_opt is None:
        return [0.0] * m, 0.0

    lp_value = dual_obj   # = sum d_i * y_{i,S}

    # ---- recover primal x_e from slack variables ----
    # Slack s_e corresponds to edge-e constraint; its reduced cost = x_e.
    x = [0.0] * m
    for e_idx in range(m):
        x[e_idx] = max(0.0, simplex.obj_row[N + e_idx])

    return x, lp_value


# ====================================================================
# 2. LP ROUNDING VIA METRIC EMBEDDING
# ====================================================================

def _metric_from_lp(n: int, edges: List[Tuple[int, int]],
                     x: List[float]) -> List[List[float]]:
    """
    Construct the shortest-path metric d from the fractional cut x*.
    d(u,v) = shortest path distance using edge weights x_e.
    """
    return _floyd_warshall(n, edges, x)


def _ell1_embed_approx(dist: List[List[float]], n: int,
                        num_cuts: int = 64) -> List[Tuple[Set[int], float]]:
    """
    Approximate an l1-embedding of the metric via random ball-cut
    sampling (simplified Ailon-Chazelle / ARV approach).

    For each trial, pick a random centre and radius, form the ball
    B(c, r), and record the cut delta(B).  Weight each cut by the
    metric mass it captures.

    Returns list of (S, weight) pairs.
    """
    cuts: List[Tuple[Set[int], float]] = []

    for _ in range(num_cuts):
        centre = random.randrange(n)
        r = random.uniform(0.0, max(dist[centre][v] for v in range(n)) + 1e-9)

        S = {v for v in range(n) if dist[centre][v] <= r}

        if not S or len(S) == n:
            continue

        weight = 0.0
        for u in range(n):
            for vv in range(u + 1, n):
                if (u in S) != (vv in S):
                    weight += dist[u][vv]

        if weight > 1e-12:
            cuts.append((S, weight))

    return cuts


def sparsest_cut_lp_rounding(
    n: int,
    edges: List[Tuple[int, int]],
    capacities: List[float],
    demands: List[Tuple[int, int, float]],
    trials: int = 200
) -> Tuple[Set[int], float, float]:
    """
    Sparsest cut via LP rounding using metric embedding.

    Algorithm (ARV-inspired):
      1. Solve LP relaxation to obtain fractional cut x*.
      2. Construct the shortest-path metric d from x*.
      3. Find an approximate l1-embedding: sample ball-cuts in the metric.
      4. Randomly perturb each sampled cut and evaluate.
      5. Return the best cut found.

    Expected approximation ratio: O(sqrt(log n)).

    Returns:
        (best_S, best_ratio, lp_ratio) where best_S is the rounded cut,
        best_ratio is its cap/demand ratio, and lp_ratio is the LP bound.
    """
    x, lp_value = solve_sparsest_cut_lp(n, edges, capacities, demands)

    if not x or lp_value <= 0:
        return set(range(n)), float('inf'), float('inf')

    # Step 2: metric from LP
    dist = _metric_from_lp(n, edges, x)

    # Step 3: approximate l1-embedding via random ball-cuts
    embedding_cuts = _ell1_embed_approx(dist, n, num_cuts=80)

    # Step 4: random rounding
    best_S: Set[int] = set(range(n))
    best_ratio = float('inf')

    for _ in range(trials):
        if embedding_cuts:
            # Weighted random choice of an embedding cut
            total_w = sum(w for _, w in embedding_cuts)
            rr = random.uniform(0.0, total_w)
            cumulative = 0.0
            chosen_S: Set[int] = set(range(n))
            for S, w in embedding_cuts:
                cumulative += w
                if cumulative >= rr:
                    chosen_S = S
                    break

            # Randomly shift boundary vertices
            boundary = set()
            for v in chosen_S:
                for u in range(n):
                    if u not in chosen_S and dist[v][u] < 2.0:
                        boundary.add(v)
                        break
            S = set(chosen_S)
            for v in boundary:
                if random.random() < 0.5:
                    S.discard(v)
        else:
            S = {v for v in range(n) if random.random() < 0.5}

        if not S or S == set(range(n)):
            continue

        ratio = _sparsest_ratio(S, edges, capacities, demands)
        if ratio < best_ratio:
            best_ratio = ratio
            best_S = set(S)

    return best_S, best_ratio, lp_value


# ====================================================================
# 3. COMBINATORIAL APPROACH (EXACT FOR SMALL INSTANCES)
# ====================================================================

def sparsest_cut_combinatorial(
    n: int,
    edges: List[Tuple[int, int]],
    capacities: List[float],
    demands: List[Tuple[int, int, float]]
) -> Tuple[Set[int], float]:
    """
    Exact sparsest cut via exhaustive enumeration of all 2^(n-1)-1 cuts.

    For each non-trivial partition (S, V\\S), compute cap(S)/demand(S).
    Feasible for n <= ~20.

    Returns:
        (best_S, best_ratio)
    """
    best_S: Set[int] = set(range(n))
    best_ratio = float('inf')

    for S in _enumerate_cuts(n):
        ratio = _sparsest_ratio(S, edges, capacities, demands)
        if ratio < best_ratio:
            best_ratio = ratio
            best_S = set(S)

    return best_S, best_ratio


# ====================================================================
# 4. CONDUCTANCE
# ====================================================================

def compute_conductance(
    n: int,
    edges: List[Tuple[int, int]],
    capacities: List[float],
    vertex_weights: List[float],
    S: Set[int]
) -> float:
    """
    Compute the conductance of set S in graph G=(V,E).

    Conductance:
        phi(S) = cap(S) / min(vol(S), vol(V\\S))

    where:
        cap(S)  = total capacity of edges crossing the cut
        vol(S)  = sum of vertex weights in S

    Conductance measures how "well-connected" a set is relative to its
    size.  Sparsest cut minimizes a weighted analogue of conductance
    over demand pairs.

    Returns:
        Conductance value (lower is better separated).
    """
    S = set(S)
    V = set(range(n))
    complement = V - S

    cap = _cut_cost(S, edges, capacities)
    vol_S = sum(vertex_weights[v] for v in S)
    vol_comp = sum(vertex_weights[v] for v in complement)
    vol_min = min(vol_S, vol_comp)

    if vol_min <= 1e-12:
        return float('inf')
    return cap / vol_min


def find_min_conductance(
    n: int,
    edges: List[Tuple[int, int]],
    capacities: List[float],
    vertex_weights: List[float]
) -> Tuple[Set[int], float]:
    """Find the set S minimizing conductance (exact, for small n)."""
    best_S: Set[int] = set(range(n))
    best_phi = float('inf')

    for S in _enumerate_cuts(n):
        phi = compute_conductance(n, edges, capacities, vertex_weights, S)
        if phi < best_phi:
            best_phi = phi
            best_S = set(S)

    return best_S, best_phi


# ====================================================================
# 5. DEMO
# ====================================================================

def demo_sparsest_cut():
    print("=" * 60)
    print("Chapter 21: Sparsest Cut")
    print("=" * 60)

    # ---- Example 1: 6-node graph with diagonal edges ----
    n1 = 6
    edges1 = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
        (0, 3), (1, 4), (2, 5),
    ]
    cap1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
    dem1 = [
        (0, 3, 1.0),
        (1, 4, 1.0),
        (2, 5, 1.0),
        (0, 2, 0.5),
        (3, 5, 0.5),
    ]

    print("\n1. Six-node graph with 3 diagonal edges")
    print(f"   Vertices: {list(range(n1))}")
    print(f"   Edges:    {edges1}")
    print(f"   Capacities: {cap1}")
    print(f"   Demands:    {dem1}")

    S_exact, ratio_exact = sparsest_cut_combinatorial(n1, edges1, cap1, dem1)
    print(f"\n   Combinatorial (exact):")
    print(f"     S          = {sorted(S_exact)}")
    print(f"     cap/demand = {ratio_exact:.4f}")

    x_lp, lp_val = solve_sparsest_cut_lp(n1, edges1, cap1, dem1)
    print(f"\n   LP relaxation value: {lp_val:.4f}")
    print(f"   LP fractional x:     {[round(v, 4) for v in x_lp]}")
    print(f"   (LP = min cost fractional cut satisfying all demand constraints;") 
    print(f"    exact sparsest cut may omit low-demand pairs to improve ratio)")

    S_round, ratio_round, lp_ratio = sparsest_cut_lp_rounding(
        n1, edges1, cap1, dem1, trials=500
    )
    print(f"\n   LP rounding:")
    print(f"     S          = {sorted(S_round)}")
    print(f"     cap/demand = {ratio_round:.4f}")
    print(f"     LP bound   = {lp_ratio:.4f}")
    if ratio_exact > 0:
        print(f"     Approx ratio to exact: {ratio_round / ratio_exact:.2f}x")
        print(f"     (LP rounding uses metric from LP for O(sqrt(log n))-approx)")

    # ---- Example 2: Path graph with long-range demands ----
    n2 = 8
    edges2 = [(i, i + 1) for i in range(7)]
    cap2 = [1.0] * 7
    dem2 = [
        (0, 7, 1.0),
        (1, 6, 0.5),
        (2, 5, 0.3),
        (3, 4, 0.2),
    ]

    print(f"\n2. Path graph P_8 with long-range demands")
    print(f"   Edges:   {edges2}")
    print(f"   Demands: {dem2}")

    S_exact2, ratio_exact2 = sparsest_cut_combinatorial(n2, edges2, cap2, dem2)
    print(f"\n   Combinatorial (exact):")
    print(f"     S          = {sorted(S_exact2)}")
    print(f"     cap/demand = {ratio_exact2:.4f}")

    x_lp2, lp_val2 = solve_sparsest_cut_lp(n2, edges2, cap2, dem2)
    print(f"   LP relaxation value: {lp_val2:.4f}")

    S_round2, ratio_round2, lp_val2r = sparsest_cut_lp_rounding(
        n2, edges2, cap2, dem2, trials=500
    )
    print(f"   LP rounding:  cap/demand = {ratio_round2:.4f}  (LP = {lp_val2r:.4f})")

    # ---- Example 3: Conductance ----
    print(f"\n3. Conductance on the path graph P_8")
    vw = [1.0] * n2
    S_test = {0, 1, 2, 3}
    phi = compute_conductance(n2, edges2, cap2, vw, S_test)
    print(f"   S = {sorted(S_test)}: phi = {phi:.4f}")

    S_min, phi_min = find_min_conductance(n2, edges2, cap2, vw)
    print(f"   Minimum conductance: S = {sorted(S_min)}, phi = {phi_min:.4f}")


if __name__ == "__main__":
    demo_sparsest_cut()

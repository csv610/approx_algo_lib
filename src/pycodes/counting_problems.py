"""
Chapter 28: Counting Problems
=============================
Vazirani Ch. 28: Counting solutions to NP problems.
- Counting DNF satisfying assignments (Karp-Luby algorithm)
- Network reliability estimation via Monte Carlo
- Near-minimum cut enumeration
"""

import math
import random
from itertools import product
from typing import List, Dict, Set, Tuple, Optional


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
Clause = List[Tuple[int, bool]]   # list of (var_index, is_positive)
Graph = Dict[int, List[int]]      # adjacency list


# ---------------------------------------------------------------------------
# 1. DNF counting  –  brute-force exact (small instances only)
# ---------------------------------------------------------------------------
def count_dnf_exact(dnf_clauses: List[Clause], n_vars: int) -> int:
    """
    Count the exact number of satisfying assignments for a DNF formula.

    Parameters
    ----------
    dnf_clauses : list of Clause
        Each clause is a list of (var_index, is_positive) tuples.
        The DNF is the disjunction of these conjunctions.
    n_vars : int
        Number of Boolean variables (indices 0 .. n_vars-1).

    Returns
    -------
    int
        Number of assignments that satisfy the DNF.

    Time complexity: O(2^n * m * k) where m = #clauses, k = max clause length.
    Only practical for n_vars <= ~20.
    """
    if n_vars > 22:
        raise ValueError(
            f"n_vars={n_vars} too large for brute-force; use Karp-Luby estimator."
        )

    count = 0
    for assignment in product((False, True), repeat=n_vars):
        if _dnf_satisfies(dnf_clauses, assignment):
            count += 1
    return count


def _dnf_satisfies(dnf_clauses: List[Clause], assignment: Tuple[bool, ...]) -> bool:
    """Check whether a given assignment satisfies the DNF formula."""
    for clause in dnf_clauses:
        if all(assignment[var] == pos for var, pos in clause):
            return True
    return False


# ---------------------------------------------------------------------------
# 2. DNF counting  –  Karp-Luby (1±ε)-approximation
# ---------------------------------------------------------------------------
def _clause_is_satisfiable(clause: Clause) -> bool:
    """Check whether a clause has no conflicting literals."""
    binding: Dict[int, bool] = {}
    for var, pos in clause:
        if var in binding:
            if binding[var] != pos:
                return False
        else:
            binding[var] = pos
    return True


def _deduplicate_clause(clause: Clause) -> Clause:
    """Remove duplicate literal bindings, keeping first occurrence."""
    seen: Dict[int, bool] = {}
    result: Clause = []
    for var, pos in clause:
        if var not in seen:
            seen[var] = pos
            result.append((var, pos))
    return result


def _sample_satisfying_assignment(
    dnf_clauses: List[Clause], n_vars: int
) -> Optional[Tuple[bool, ...]]:
    """
    Sample a random assignment that satisfies exactly one uniformly chosen
    clause.  If the formula is unsatisfiable, return None.

    For clause C_i of length k_i (after deduplication):
    - Assign the k_i literal variables to make C_i true.
    - Assign the remaining n_vars - k_i variables uniformly at random.
    - Probability of this assignment: (1/2)^(n_vars - k_i).

    Unsatisfiable clauses (with conflicting literals) are skipped.
    """
    # Filter to satisfiable clauses and pick one uniformly
    sat_clauses = [c for c in dnf_clauses if _clause_is_satisfiable(c)]
    if not sat_clauses:
        return None

    clause = random.choice(sat_clauses)
    clause = _deduplicate_clause(clause)

    assignment = [False] * n_vars
    clause_vars = set()

    for var, pos in clause:
        assignment[var] = pos
        clause_vars.add(var)

    for v in range(n_vars):
        if v not in clause_vars:
            assignment[v] = random.random() < 0.5

    return tuple(assignment)


def _assignment_prob(dnf_clauses: List[Clause], assignment: Tuple[bool, ...]) -> float:
    """
    Compute the probability that _sample_satisfying_assignment produces the
    given assignment (which must satisfy the DNF).

    Only satisfiable clauses (no conflicting literals) participate in sampling.
    For each satisfiable clause C_i satisfied by the assignment, let k_i be
    the number of distinct variables in C_i.  The assignment is produced by
    choosing clause i with probability 1/m' (where m' = #satisfiable clauses)
    and then randomising n - k_i free variables, so the probability
    contribution from clause i is:

        (1/m') * (1/2)^{n - k_i}    if assignment agrees with C_i's literals

    We sum over all satisfiable clauses whose literals the assignment satisfies.
    """
    sat_clauses = [c for c in dnf_clauses if _clause_is_satisfiable(c)]
    m = len(sat_clauses)
    if m == 0:
        return 0.0

    n_vars = len(assignment)
    total_prob = 0.0

    for clause in sat_clauses:
        deduped = _deduplicate_clause(clause)
        if all(assignment[var] == pos for var, pos in deduped):
            k_i = len(deduped)
            total_prob += (1.0 / m) * (0.5 ** (n_vars - k_i))

    return total_prob


def count_dnf_karp_luby(
    dnf_clauses: List[Clause],
    n_vars: int,
    epsilon: float = 0.1,
    trials: int = 1000,
) -> Tuple[float, float, float]:
    """
    Karp-Luby estimator for the number of satisfying assignments of a DNF.

    Returns a (1 ± ε)-approximation with high probability.

    Algorithm
    ---------
    1. Repeatedly sample a random satisfying assignment by choosing a clause
       uniformly, fixing its literals, and randomising the rest.
    2. For each sample x_j, compute f_j = 1 / Pr[x_j] where Pr[x_j] is
       the probability that x_j would be produced by the sampling procedure.
    3. The estimator is  T = (1/trials) * sum_j f_j * I[x_j satisfies DNF].
       Since every sample satisfies the DNF by construction, T is an
       unbiased estimator of |SAT(DNF)|.

    Parameters
    ----------
    dnf_clauses : list of Clause
    n_vars : int
    epsilon : float
        Desired relative error.  With high probability the estimate
        lies in [(1-ε) * exact, (1+ε) * exact].
    trials : int
        Number of Monte Carlo samples.  O(m / ε²) suffices for
        a (1±ε)-approximation w.h.p., where m = #clauses.

    Returns
    -------
    (estimate, lower_bound, upper_bound)
        The point estimate and the confidence interval.
    """
    if not dnf_clauses:
        return 0.0, 0.0, 0.0

    m = len(dnf_clauses)
    total = 0.0
    total_sq = 0.0

    for _ in range(trials):
        sample = _sample_satisfying_assignment(dnf_clauses, n_vars)
        if sample is None:
            f = 0.0
        else:
            prob = _assignment_prob(dnf_clauses, sample)
            f = 1.0 / prob if prob > 0 else 0.0
        total += f
        total_sq += f * f

    estimate = total / trials

    # Variance estimate
    mean = total / trials
    variance = total_sq / trials - mean * mean
    std_err = math.sqrt(max(variance, 0) / trials)

    # 95 % confidence interval
    lower = estimate - 1.96 * std_err
    upper = estimate + 1.96 * std_err

    return estimate, max(lower, 0.0), upper


# ---------------------------------------------------------------------------
# 3. Network reliability  –  brute-force exact
# ---------------------------------------------------------------------------
def network_reliability_exact(
    graph: Graph, failure_probs: Dict[Tuple[int, int], float]
) -> float:
    """
    Compute the exact probability that the graph remains connected when each
    edge independently fails with the given probability.

    Parameters
    ----------
    graph : Graph
        Adjacency list.  Undirected: (u,v) and (v,u) both present.
    failure_probs : dict mapping (min(u,v), max(u,v)) -> probability of failure.

    Returns
    -------
    float
        Probability that the surviving subgraph is connected.

    Time complexity: O(2^E * poly(n)).  Only practical for E <= ~25.
    """
    # Collect all edges (canonical form)
    edges = set()
    for u in graph:
        for v in graph[u]:
            edge = (min(u, v), max(u, v))
            edges.add(edge)
    edges = sorted(edges)
    n_edges = len(edges)

    if n_edges > 25:
        raise ValueError(
            f"n_edges={n_edges} too large for brute-force; use Monte Carlo estimator."
        )

    nodes = set(graph.keys())
    reliability = 0.0

    # Enumerate all subsets of surviving edges
    for mask in range(1 << n_edges):
        prob = 1.0
        surviving = {n: [] for n in nodes}

        for i, (u, v) in enumerate(edges):
            if mask & (1 << i):
                # Edge survives
                prob *= 1.0 - failure_probs.get((u, v), 0.0)
                surviving[u].append(v)
                surviving[v].append(u)
            else:
                # Edge fails
                prob *= failure_probs.get((u, v), 0.0)

        if _is_connected(surviving, nodes):
            reliability += prob

    return reliability


def _is_connected(adj: Graph, nodes: Set[int]) -> bool:
    """Check if the graph (given as adjacency list) is connected via BFS."""
    if not nodes:
        return True
    start = next(iter(nodes))
    visited = {start}
    queue = [start]
    while queue:
        u = queue.pop()
        for v in adj.get(u, []):
            if v not in visited:
                visited.add(v)
                queue.append(v)
    return visited == nodes


# ---------------------------------------------------------------------------
# 4. Network reliability  –  Monte Carlo estimator
# ---------------------------------------------------------------------------
def network_reliability_monte_carlo(
    graph: Graph,
    failure_probs: Dict[Tuple[int, int], float],
    trials: int = 10000,
) -> Tuple[float, float, float]:
    """
    Monte Carlo estimator for network reliability.

    Each trial independently deletes every edge with its failure probability,
    then checks connectivity of the surviving graph.

    Returns
    -------
    (estimate, lower_bound, upper_bound)
        The point estimate and the 95 % confidence interval.
    """
    # Collect edges
    edges = set()
    for u in graph:
        for v in graph[u]:
            edges.add((min(u, v), max(u, v)))
    edges = sorted(edges)

    nodes = set(graph.keys())
    successes = 0

    for _ in range(trials):
        surviving = {n: [] for n in nodes}
        for u, v in edges:
            q = failure_probs.get((u, v), 0.0)
            if random.random() >= q:
                surviving[u].append(v)
                surviving[v].append(u)
        if _is_connected(surviving, nodes):
            successes += 1

    p_hat = successes / trials
    std_err = math.sqrt(p_hat * (1.0 - p_hat) / trials)

    return p_hat, max(p_hat - 1.96 * std_err, 0.0), min(p_hat + 1.96 * std_err, 1.0)


# ---------------------------------------------------------------------------
# 5. Near-minimum cut enumeration
# ---------------------------------------------------------------------------
def _global_min_cut_karger(graph: Graph, trials: int = 100) -> int:
    """
    Estimate the global minimum cut size using Karger's contraction algorithm.
    Returns the best (smallest) cut found across `trials` independent runs.
    """
    nodes = list(graph.keys())
    n = len(nodes)
    best = float('inf')

    for _ in range(trials):
        # Build edge list with multiplicity
        adj = {u: list(neighbors) for u, neighbors in graph.items()}
        num_nodes = n

        while num_nodes > 2:
            # Pick a random edge
            u = random.choice(list(adj.keys()))
            if not adj[u]:
                break
            v = random.choice(adj[u])

            # Contract: merge v into u
            adj[u] = [x for x in adj[u] if x != v] + [x for x in adj[v] if x != u]
            for w in adj[v]:
                adj[w] = [u if x == v else x for x in adj[w]]
            del adj[v]
            num_nodes -= 1

        # Count remaining edges incident to any surviving node
        remaining_edges = sum(len(neighbors) for neighbors in adj.values()) // 2
        best = min(best, remaining_edges)

    return best


def count_near_min_cuts(
    graph: Graph,
    threshold_factor: float = 1.01,
    max_cuts: int = 1000,
) -> Tuple[int, int, float]:
    """
    Count the number of near-minimum cuts in the graph.

    A cut (S, V\\S) is a *λ(n/λ)-near-minimum cut* if its size is at most
    λ · threshold_factor, where λ is the global minimum cut size.

    Uses Karger's result: the number of (λ·c)-near-min cuts is at most n^{2c}.
    We enumerate cuts by repeatedly running Karger's contraction and collecting
    distinct cuts.

    Parameters
    ----------
    graph : Graph
        Adjacency list (unweighted).  Self-loops and parallel edges ignored.
    threshold_factor : float
        Multiplier on the minimum cut.  A cut is "near" if its size ≤
        λ * threshold_factor.  Default 1.01 means within 1 % of min cut.
    max_cuts : int
        Maximum number of distinct near-min cuts to enumerate.

    Returns
    -------
    (min_cut_size, num_near_cuts, upper_bound)
        The minimum cut size, the number of distinct near-min cuts found,
        and a theoretical upper bound n^{2c} from Karger's analysis.
    """
    nodes = set(graph.keys())
    n = len(nodes)

    if n <= 1:
        return 0, 0, 0

    # Estimate min cut
    lam = _global_min_cut_karger(graph, trials=50)
    if lam == 0:
        # Disconnected graph
        return 0, float('inf'), float('inf')

    threshold = int(math.ceil(lam * threshold_factor))
    c = threshold_factor

    # Enumerate distinct near-min cuts via repeated Karger runs
    found_cuts: Set[frozenset] = set()

    for _ in range(min(max_cuts * 10, 5000)):
        if len(found_cuts) >= max_cuts:
            break
        cut = _karger_single_cut(graph)
        if cut is not None and len(cut) <= threshold and len(cut) > 0:
            found_cuts.add(frozenset(cut))

    upper_bound = n ** (2 * c) if c < n else float('inf')

    return lam, len(found_cuts), upper_bound


def _karger_single_cut(graph: Graph) -> Optional[Set[Tuple[int, int]]]:
    """
    Run a single Karger contraction and return the cut edges.
    Returns a set of (min, max) edge tuples forming the cut.
    """
    nodes = list(graph.keys())
    n = len(nodes)
    if n <= 1:
        return None

    # Build adjacency with canonical edges
    adj: Dict[int, List[int]] = {u: [] for u in graph}
    for u in graph:
        for v in graph[u]:
            adj[u].append(v)

    parent = {u: u for u in graph}
    rank = {u: 0 for u in graph}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx == ry:
            return
        if rank[rx] < rank[ry]:
            rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]:
            rank[rx] += 1

    alive = set(graph.keys())

    while len(alive) > 2:
        # Pick random edge
        u = random.choice(list(alive))
        if not adj[u]:
            break
        v = random.choice(adj[u])
        if find(u) == find(v):
            # Self-loop – retry by picking a different neighbour
            adj[u] = [x for x in adj[u] if find(x) != find(u)]
            if not adj[u]:
                break
            v = random.choice(adj[u])
            if find(u) == find(v):
                continue

        # Contract v into u
        union(u, v)
        alive.discard(v)
        # Merge adjacency
        adj[u] = [x for x in adj[u] if find(x) != find(u)] + \
                  [x for x in adj[v] if find(x) != find(u)]

    # Collect cut edges between the two surviving super-nodes
    if len(alive) != 2:
        return None

    a, b = list(alive)
    ra, rb = find(a), find(b)
    cut_edges: Set[Tuple[int, int]] = set()
    for u in graph:
        for v in graph[u]:
            if find(u) == ra and find(v) == rb:
                cut_edges.add((min(u, v), max(u, v)))
            elif find(u) == rb and find(v) == ra:
                cut_edges.add((min(u, v), max(u, v)))

    return cut_edges


# ---------------------------------------------------------------------------
# 6. Demo
# ---------------------------------------------------------------------------
def demo_counting_problems() -> None:
    """Demonstrate counting problem algorithms."""
    print("=" * 65)
    print("Chapter 28: Counting Problems")
    print("=" * 65)

    # --- DNF counting ---------------------------------------------------
    print("\n--- DNF Counting (Karp-Luby) ---")
    # Formula: (x0 AND x1) OR (NOT x1 AND x2) OR (x0 AND NOT x2)
    dnf = [
        [(0, True), (1, True)],        # x0 ∧ x1
        [(1, False), (2, True)],        # ¬x1 ∧ x2
        [(0, True), (2, False)],        # x0 ∧ ¬x2
    ]
    n_vars = 3

    exact = count_dnf_exact(dnf, n_vars)
    print(f"  DNF: (x0∧x1) ∨ (¬x1∧x2) ∨ (x0∧¬x2)")
    print(f"  Variables: {n_vars}")
    print(f"  Exact count: {exact} out of {2**n_vars} assignments")

    est, lo, hi = count_dnf_karp_luby(dnf, n_vars, epsilon=0.1, trials=2000)
    print(f"  Karp-Luby estimate: {est:.2f}  [{lo:.2f}, {hi:.2f}]")

    # Larger instance
    print("\n  Larger instance (8 vars, 12 clauses):")
    random.seed(42)
    big_dnf = []
    for _ in range(12):
        k = random.randint(2, 4)
        chosen_vars = random.sample(range(8), k)
        clause = [(v, random.random() < 0.5) for v in chosen_vars]
        big_dnf.append(clause)
    big_exact = count_dnf_exact(big_dnf, 8)
    big_est, big_lo, big_hi = count_dnf_karp_luby(big_dnf, 8, epsilon=0.05, trials=5000)
    print(f"  Exact: {big_exact} / {256}")
    print(f"  Estimate: {big_est:.2f}  [{big_lo:.2f}, {big_hi:.2f}]")
    relative_error = abs(big_est - big_exact) / big_exact if big_exact > 0 else 0
    print(f"  Relative error: {relative_error:.4f}")

    # --- Network reliability --------------------------------------------
    print("\n--- Network Reliability ---")
    #  1 -- 2 -- 3
    #  |         |
    #  4 ------- 5
    g: Graph = {
        1: [2, 4],
        2: [1, 3],
        3: [2, 5],
        4: [1, 5],
        5: [3, 4],
    }
    # Each edge fails with probability 0.1
    edges_list = [(1, 2), (1, 4), (2, 3), (3, 5), (4, 5)]
    fail_probs = {e: 0.1 for e in edges_list}

    exact_rel = network_reliability_exact(g, fail_probs)
    print(f"  Graph: 5-node cycle-like (5 edges)")
    print(f"  Edge failure probability: 0.1")
    print(f"  Exact reliability: {exact_rel:.6f}")

    mc_est, mc_lo, mc_hi = network_reliability_monte_carlo(g, fail_probs, trials=20000)
    print(f"  Monte Carlo estimate: {mc_est:.6f}  [{mc_lo:.6f}, {mc_hi:.6f}]")
    print(f"  Absolute error: {abs(mc_est - exact_rel):.6f}")

    # --- Near-minimum cuts ----------------------------------------------
    print("\n--- Near-Minimum Cuts ---")
    # Petersen-like small graph
    petersen: Graph = {
        0: [1, 4, 5],
        1: [0, 2, 6],
        2: [1, 3, 7],
        3: [2, 4, 8],
        4: [3, 0, 9],
        5: [0, 7, 8],
        6: [1, 8, 9],
        7: [2, 5, 9],
        8: [3, 5, 6],
        9: [4, 6, 7],
    }
    lam, n_near, ub = count_near_min_cuts(petersen, threshold_factor=1.5)
    print(f"  Petersen graph (10 nodes, 15 edges)")
    print(f"  Minimum cut size λ = {lam}")
    print(f"  Near-min cuts found (≤{1.5:.1f}λ): {n_near}")
    print(f"  Theoretical upper bound n^{{2c}}: {ub:.0f}")

    print("\n" + "=" * 65)


if __name__ == "__main__":
    demo_counting_problems()

"""
Chapter 29: Hardness of Approximation (PCP Theorem)
====================================================
Vazirani Ch. 29: PCP theorem and hardness of approximation results.
Implements:
1. PCP verifier simulation with O(log n) random bits and O(1) queries.
2. Random 3-SAT instance generation for MAX-3SAT hardness demonstration.
3. Brute-force exact MAX-3SAT solver for small instances.
4. Greedy MAX-3SAT heuristic.
5. Hardness gap demonstration showing the 7/8 bound for random 3-SAT.
"""

import math
import random
from typing import List, Tuple, Set, Optional


# ---------------------------------------------------------------------------
# PCP Verifier Simulation
# ---------------------------------------------------------------------------

def _binary_string_to_int(bits: List[int]) -> int:
    """Convert a list of 0/1 bits to an integer (big-endian)."""
    val = 0
    for b in bits:
        val = (val << 1) | b
    return val


def verify_pcp(
    proof_bits: List[int],
    n_random_bits: int = 4,
    n_queries: int = 3,
    verifier_seed: Optional[int] = None,
) -> Tuple[bool, List[Tuple[int, int]], float]:
    """
    Simulate a PCP verifier for a decision problem.

    The verifier:
    1. Picks random bits r (n_random_bits of them).
    2. Uses r to compute q query indices into the proof.
    3. Accepts iff the queried proof bits satisfy a fixed predicate.

    For MAX-3SAT the predicate is a random 3-CNF clause on the q queried
    positions (interpreted as a 3-bit index when q=3).

    Returns:
        accepted:  Whether the verifier accepted.
        queries:   List of (index, bit_value) pairs read from the proof.
        prob_accept: Estimated acceptance probability (over many random seeds).
    """
    rng = random.Random(verifier_seed)

    if len(proof_bits) == 0:
        return False, [], 0.0

    n = len(proof_bits)

    # --- Run the verifier once with the given (or fresh) randomness -------
    r_bits = [rng.randint(0, 1) for _ in range(n_random_bits)]
    r_int = _binary_string_to_int(r_bits)

    # Deterministic query selection: hash r_int to q indices
    query_indices = []
    for q in range(n_queries):
        idx = (r_int * (q + 7) + q * 13) % n
        query_indices.append(idx)

    queries = [(idx, proof_bits[idx]) for idx in query_indices]

    # --- Predicate: majority of queried bits must be 1 (simple PCP check) --
    ones = sum(bit for _, bit in queries)
    accepted = ones >= math.ceil(n_queries / 2)

    # --- Estimate acceptance probability over many random seeds -----------
    total_trials = 256
    accept_count = 0
    for trial in range(total_trials):
        trial_rng = random.Random(verifier_seed if verifier_seed is not None else trial)
        t_rbits = [trial_rng.randint(0, 1) for _ in range(n_random_bits)]
        t_rint = _binary_string_to_int(t_rbits)
        t_indices = []
        for q in range(n_queries):
            idx = (t_rint * (q + 7) + q * 13) % n
            t_indices.append(idx)
        t_ones = sum(proof_bits[i] for i in t_indices)
        if t_ones >= math.ceil(n_queries / 2):
            accept_count += 1

    return accepted, queries, accept_count / total_trials


# ---------------------------------------------------------------------------
# MAX-3SAT Instance Generation and Solvers
# ---------------------------------------------------------------------------

def max3sat_random_instance(
    n_vars: int,
    n_clauses: int,
    seed: Optional[int] = None,
) -> List[Tuple[Tuple[int, bool], Tuple[int, bool], Tuple[int, bool]]]:
    """
    Generate a random 3-SAT instance.

    Each clause has exactly 3 literals. A literal is (var_index, is_positive).
    var_index is 0-based.

    Returns:
        clauses: List of 3-tuples of literals.
    """
    rng = random.Random(seed)
    clauses = []
    for _ in range(n_clauses):
        lits = []
        for _ in range(3):
            var = rng.randint(0, n_vars - 1)
            pos = rng.choice([True, False])
            lits.append((var, pos))
        clauses.append(tuple(lits))
    return clauses


def _evaluate_clause(
    clause: Tuple[Tuple[int, bool], Tuple[int, bool], Tuple[int, bool]],
    assignment: List[bool],
) -> bool:
    """Check if a single 3-SAT clause is satisfied by the assignment."""
    for var, is_pos in clause:
        if assignment[var] == is_pos:
            return True
    return False


def max3sat_bruteforce(
    clauses: List[Tuple[Tuple[int, bool], Tuple[int, bool], Tuple[int, bool]]],
    n_vars: Optional[int] = None,
) -> Tuple[int, List[bool], float]:
    """
    Exact MAX-3SAT solver via brute-force enumeration.

    Tries all 2^n_vars assignments and returns the best.
    Only practical for small n_vars (<= ~20).

    Returns:
        max_sat:       Number of satisfied clauses.
        best_assign:   Boolean assignment list (0-indexed).
        ratio:         max_sat / len(clauses).
    """
    if n_vars is None:
        n_vars = max(var for clause in clauses for var, _ in clause) + 1

    best_sat = -1
    best_assign: List[bool] = [False] * n_vars

    for mask in range(1 << n_vars):
        assign = [(mask >> i) & 1 == 1 for i in range(n_vars)]
        sat = sum(1 for c in clauses if _evaluate_clause(c, assign))
        if sat > best_sat:
            best_sat = sat
            best_assign = assign

    ratio = best_sat / len(clauses) if clauses else 0.0
    return best_sat, best_assign, ratio


def max3sat_greedy(
    clauses: List[Tuple[Tuple[int, bool], Tuple[int, bool], Tuple[int, bool]]],
    n_vars: Optional[int] = None,
    n_flips: int = 100,
    seed: Optional[int] = None,
) -> Tuple[int, List[bool], float]:
    """
    Greedy / local-search MAX-3SAT heuristic.

    Starts from a random assignment, then repeatedly flips the variable
    that gives the largest improvement (or the least deterioration if no
    improvement exists).  Stops after n_flips with no improvement.

    Returns:
        sat_count:    Number of satisfied clauses.
        assignment:   Boolean assignment list (0-indexed).
        ratio:        sat_count / len(clauses).
    """
    if not clauses:
        return 0, [], 1.0

    if n_vars is None:
        n_vars = max(var for clause in clauses for var, _ in clause) + 1

    rng = random.Random(seed)
    assign = [rng.choice([True, False]) for _ in range(n_vars)]

    def count_sat(a: List[bool]) -> int:
        return sum(1 for c in clauses if _evaluate_clause(c, a))

    best_sat = count_sat(assign)
    no_improve = 0

    while no_improve < n_flips:
        best_flip = -1
        best_after = best_sat

        for v in range(n_vars):
            assign[v] = not assign[v]
            s = count_sat(assign)
            if s > best_after:
                best_after = s
                best_flip = v
            assign[v] = not assign[v]

        if best_flip >= 0:
            assign[best_flip] = not assign[best_flip]
            best_sat = best_after
            no_improve = 0
        else:
            no_improve += 1

    ratio = best_sat / len(clauses) if clauses else 0.0
    return best_sat, assign, ratio


# ---------------------------------------------------------------------------
# Hardness Gap Demonstration
# ---------------------------------------------------------------------------

def demonstrate_hardness_gap(
    n_vars: int = 8,
    n_clauses: int = 40,
    n_trials: int = 5,
    seed: int = 42,
) -> List[Tuple[int, int, float, float, float]]:
    """
    Show the hardness gap for MAX-3SAT.

    For random 3-CNF with clauses of 3 literals, the expected fraction of
    clauses satisfiable by a random assignment is 7/8.  The PCP theorem
    implies that approximating MAX-3SAT within factor 7/8 + epsilon is
    NP-hard (for any epsilon > 0).

    This function:
    1. Generates random 3-SAT instances.
    2. Solves each exactly (brute-force).
    3. Solves each with the greedy heuristic.
    4. Compares both to the 7/8 random-assignment baseline.

    Returns:
        List of (n_vars, n_clauses, exact_ratio, greedy_ratio, random_ratio)
        tuples for each trial.
    """
    results = []
    for i in range(n_trials):
        clauses = max3sat_random_instance(n_vars, n_clauses, seed=seed + i)
        exact_sat, _, exact_ratio = max3sat_bruteforce(clauses, n_vars)
        greedy_sat, _, greedy_ratio = max3sat_greedy(clauses, n_vars, seed=seed + i)

        # Random assignment expected ratio
        random_ratio = 7.0 / 8.0

        results.append((n_vars, n_clauses, exact_ratio, greedy_ratio, random_ratio))

    return results


def verify_78_baseline(n_vars: int = 10, n_clauses: int = 200, n_trials: int = 50, seed: int = 0) -> float:
    """
    Empirically verify that a random assignment satisfies 7/8 of random
    3-SAT clauses on average.

    Returns:
        Average fraction of clauses satisfied.
    """
    rng = random.Random(seed)
    total_ratios = []

    for _ in range(n_trials):
        clauses = max3sat_random_instance(n_vars, n_clauses, seed=rng.randint(0, 10**9))
        assign = [rng.choice([True, False]) for _ in range(n_vars)]
        sat = sum(1 for c in clauses if _evaluate_clause(c, assign))
        total_ratios.append(sat / n_clauses)

    return sum(total_ratios) / len(total_ratios)


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo_hardness_approx() -> None:
    print("=" * 65)
    print("Chapter 29: Hardness of Approximation — PCP Theorem")
    print("=" * 65)

    # --- 1. PCP Verifier Demo ---
    print("\n1. PCP Verifier Simulation")
    print("   A PCP verifier reads O(1) bits from a proof string using O(log n)")
    print("   random bits and decides to accept/reject.\n")

    n_proof = 16
    proof_all_ones = [1] * n_proof
    proof_all_zeros = [0] * n_proof
    proof_mixed = [random.choice([0, 1]) for _ in range(n_proof)]

    for label, proof in [("all-ones proof", proof_all_ones),
                         ("all-zeros proof", proof_all_zeros),
                         ("random proof", proof_mixed)]:
        acc, queries, p_accept = verify_pcp(proof, n_random_bits=4, n_queries=3, verifier_seed=7)
        print(f"   {label:20s} -> accepted={acc}, "
              f"queries={queries}, P(accept)={p_accept:.3f}")

    # --- 2. 7/8 Baseline Verification ---
    print("\n2. Random Assignment Baseline for 3-SAT (expected 7/8 = 0.8750)")
    avg = verify_78_baseline(n_vars=10, n_clauses=200, n_trials=50, seed=42)
    print(f"   Empirical average fraction satisfied: {avg:.4f}")
    print("   (This is the baseline: random 3-SAT has a 7/8 satisfiable fraction.)")

    # --- 3. MAX-3SAT Exact vs Greedy vs 7/8 Baseline ---
    print("\n3. Hardness Gap Demonstration (MAX-3SAT)")
    print(f"   {'n_vars':>6s} {'n_clauses':>9s} {'exact':>8s} {'greedy':>8s} {'7/8 baseline':>12s}")
    print("   " + "-" * 50)

    for nv, nc in [(5, 20), (6, 30), (7, 40), (8, 50)]:
        results = demonstrate_hardness_gap(n_vars=nv, n_clauses=nc, n_trials=3, seed=10)
        for _, _, er, gr, rr in results:
            print(f"   {nv:6d} {nc:9d} {er:8.4f} {gr:8.4f} {rr:12.4f}")

    # --- 4. Hardness Implications ---
    print("\n4. Hardness of Approximation Summary (Vazirani Ch. 29)")
    print("   -----------------------------------------------------------")
    print("   Problem            | Hardness result (unless P=NP)")
    print("   -----------------------------------------------------------")
    print("   MAX-3SAT            | Cannot approx within 7/8 + eps")
    print("   Vertex Cover        | Cannot approx within 2 - eps (UGC)")
    print("   Clique              | No constant-factor approx")
    print("   Set Cover           | Cannot approx within (1-eps) ln n")
    print("   -----------------------------------------------------------")
    print("   PCP Theorem: NP = PCP(O(log n), O(1))")
    print("   Every NP problem has a proof checkable with O(log n) random")
    print("   bits and O(1) query bits.")

    # --- 5. Small Instance: Exact Optimal vs Greedy ---
    print("\n5. Small Instance Detail (5 vars, 20 clauses)")
    clauses = max3sat_random_instance(5, 20, seed=99)
    exact_sat, exact_assign, exact_ratio = max3sat_bruteforce(clauses, 5)
    greedy_sat, greedy_assign, greedy_ratio = max3sat_greedy(clauses, 5, seed=99)
    print(f"   Clauses: {len(clauses)}")
    print(f"   Exact:   {exact_sat}/{len(clauses)} satisfied ({exact_ratio:.4f})")
    print(f"            assignment = {exact_assign}")
    print(f"   Greedy:  {greedy_sat}/{len(clauses)} satisfied ({greedy_ratio:.4f})")
    print(f"            assignment = {greedy_assign}")
    print(f"   7/8 baseline: {7/8:.4f}")


if __name__ == "__main__":
    demo_hardness_approx()

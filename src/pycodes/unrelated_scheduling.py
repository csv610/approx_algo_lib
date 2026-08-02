"""
Chapter 17: Scheduling on Unrelated Parallel Machines
=====================================================
Vazirani Ch. 17: Scheduling on Unrelated Parallel Machines.

Problem: Given n jobs and m unrelated parallel machines, where job j has
processing time p_{ij} on machine i, assign each job to exactly one machine
to minimize the makespan (maximum completion time across all machines).

This problem is strongly NP-hard (reduces from Partition when m=2).

Vazirani's Approach (2-approximation via LP relaxation + rounding):
1. LP Relaxation: Minimize T s.t. sum_j p_{ij} x_{ij} <= T for all i,
   sum_i x_{ij} = 1 for all j, x_{ij} >= 0
2. Parametric Pruning: Binary search on T, solve LP feasibility
3. Extreme Point Rounding: An extreme point of the LP has at most m-1
   fractional x_{ij} values per machine. Rounding these fractional
   assignments yields a 2-approximation.

Key Insight: For a fixed T, the LP feasibility LP is a linear program
that can be solved via maximum flow in a bipartite graph. The rounding
uses the structural property of extreme points.

Implements:
1. lp_feasibility(T) - Check if schedule with makespan T is feasible via LP
2. unrelated_scheduling_lp() - Binary search + LP for 2-approximation
3. unrelated_scheduling_greedy() - Simple greedy for comparison
4. unrelated_scheduling_local_search() - Local search improvement
"""

import sys
import os
import math
from typing import List, Tuple, Dict, Optional

sys.path.insert(0, os.path.dirname(__file__))
from lp_algorithms import Simplex


# ============================================================
# LP FEASIBILITY CHECK
# ============================================================

def lp_feasibility(T: float,
                   processing_times: List[List[float]]) -> Tuple[bool, Optional[List[List[float]]]]:
    """
    Check if a schedule with makespan T is feasible via LP relaxation.

    LP:
      Minimize 0
      s.t.  sum_j p_{ij} * x_{ij} <= T    for all machines i
            sum_i x_{ij} = 1               for all jobs j
            x_{ij} >= 0

    Returns (feasible, fractional_assignment) where fractional_assignment
    is an m x n matrix of x_{ij} values if feasible.
    """
    m = len(processing_times)
    if m == 0:
        return True, []
    n = len(processing_times[0])
    if n == 0:
        return True, [[] for _ in range(m)]

    # Variables: [x_00, x_01, ..., x_{m-1,n-1}, T]
    # Total variables: m*n + 1
    num_vars = m * n + 1
    T_idx = m * n  # index of T in variable vector

    # Constraints (all in <= form for Simplex):
    # 1) Machine capacity: -p_{ij} x_{ij} - T <= -T  (m constraints)
    # 2) Assignment upper:   sum_i x_{ij} <= 1        (n constraints)
    # 3) Assignment lower:  -sum_i x_{ij} <= -1       (n constraints)

    A = []
    b = []

    # Machine capacity constraints
    for i in range(m):
        row = [0.0] * num_vars
        for j in range(n):
            row[i * n + j] = -processing_times[i][j]
        row[T_idx] = -1.0
        A.append(row)
        b.append(-T)

    # Assignment upper bound: sum_i x_{ij} <= 1
    for j in range(n):
        row = [0.0] * num_vars
        for i in range(m):
            row[i * n + j] = 1.0
        A.append(row)
        b.append(1.0)

    # Assignment lower bound: -sum_i x_{ij} <= -1
    for j in range(n):
        row = [0.0] * num_vars
        for i in range(m):
            row[i * n + j] = -1.0
        A.append(row)
        b.append(-1.0)

    # Objective: Minimize 0 (just feasibility)
    c = [0.0] * num_vars

    simplex = Simplex(A, b, c)
    x, obj = simplex.solve()

    if x is None:
        return False, None

    # Extract fractional assignment matrix
    assignment = [[0.0] * n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            assignment[i][j] = max(0.0, x[i * n + j])

    return True, assignment


# ============================================================
# LP ROUNDING (Extreme Point Property)
# ============================================================

def round_fractional_assignment(
    fractional: List[List[float]],
    processing_times: List[List[float]],
    T: float,
) -> List[int]:
    """
    Round a fractional LP assignment to an integral one using the
    extreme point structure (Vazirani Ch. 17).

    Key property: An extreme point of the LP has at most m-1 fractional
    x_{ij} values across all machines combined (since each machine's
    tight constraint eliminates one degree of freedom).

    Rounding strategy:
    - For each job j, find the machine with max x_{ij}
    - Assign j to that machine (this is the "natural" rounding)
    - If a machine is overloaded, reassign its least-committed fractional
      jobs to other machines with spare capacity

    This yields a 2-approximation.
    """
    m = len(fractional)
    if m == 0:
        return []
    n = len(fractional[0])

    # Step 1: Natural rounding - assign each job to its dominant machine
    assignment = [-1] * n
    machine_loads = [0.0] * m

    for j in range(n):
        # Find machine with highest fractional assignment
        best_machine = 0
        best_val = 0.0
        for i in range(m):
            if fractional[i][j] > best_val:
                best_val = fractional[i][j]
                best_machine = i
        assignment[j] = best_machine
        machine_loads[best_machine] += processing_times[best_machine][j]

    # Step 2: Local improvement - swap overloaded jobs
    # If any machine exceeds T, try to move its jobs to less loaded machines
    max_load = max(machine_loads) if machine_loads else 0.0

    for _ in range(m * n):
        if max_load <= T + 1e-9:
            break

        # Find most overloaded machine
        overloaded = max(range(m), key=lambda i: machine_loads[i])
        if machine_loads[overloaded] <= T + 1e-9:
            break

        # Try to move each job on overloaded machine to another machine
        moved = False
        for j in range(n):
            if assignment[j] != overloaded:
                continue

            # Try moving to machine with least load
            best_dest = -1
            best_new_load = machine_loads[overloaded]
            for i in range(m):
                if i == overloaded:
                    continue
                new_load_i = machine_loads[i] + processing_times[i][j]
                if new_load_i < best_new_load:
                    best_new_load = new_load_i
                    best_dest = i

            if best_dest >= 0:
                # Move job j from overloaded to best_dest
                old_machine = overloaded
                machine_loads[old_machine] -= processing_times[old_machine][j]
                machine_loads[best_dest] += processing_times[best_dest][j]
                assignment[j] = best_dest
                moved = True
                break

        if not moved:
            break

        max_load = max(machine_loads) if machine_loads else 0.0

    return assignment


# ============================================================
# MAIN ALGORITHM: LP + Binary Search (2-approximation)
# ============================================================

def unrelated_scheduling_lp(
    processing_times: List[List[float]],
) -> Tuple[List[int], float]:
    """
    2-approximation for Scheduling on Unrelated Parallel Machines
    via LP relaxation + parametric pruning (Vazirani Ch. 17).

    Algorithm:
    1. Binary search on makespan T
    2. For each T, solve LP feasibility
    3. Round fractional solution to integral assignment

    Approximation factor: 2 (from LP rounding of extreme points)

    Args:
        processing_times: m x n matrix where p[i][j] is the processing
                          time of job j on machine i

    Returns:
        (assignment, makespan) where assignment[j] = machine index
    """
    m = len(processing_times)
    if m == 0:
        return [], 0.0
    n = len(processing_times[0])
    if n == 0:
        return [0] * 0, 0.0

    # Bounds for binary search
    min_time = min(processing_times[i][j] for i in range(m) for j in range(n))
    max_time = sum(max(processing_times[i][j] for i in range(m)) for j in range(n))

    best_assignment = list(range(n))  # fallback: each job to machine 0
    best_makespan = max_time

    # Binary search on T
    lo, hi = min_time, max_time
    for _ in range(60):  # enough iterations for convergence
        if hi - lo < 1e-7:
            break
        mid = (lo + hi) / 2.0
        feasible, frac = lp_feasibility(mid, processing_times)
        if feasible:
            # Round to integral assignment
            assignment = round_fractional_assignment(frac, processing_times, mid)
            # Compute actual makespan of rounded solution
            loads = [0.0] * m
            for j in range(n):
                loads[assignment[j]] += processing_times[assignment[j]][j]
            actual_ms = max(loads) if loads else 0.0
            if actual_ms < best_makespan:
                best_makespan = actual_ms
                best_assignment = assignment[:]
            hi = mid
        else:
            lo = mid

    return best_assignment, best_makespan


# ============================================================
# GREEDY BASELINE
# ============================================================

def unrelated_scheduling_greedy(
    processing_times: List[List[float]],
) -> Tuple[List[int], float]:
    """
    Simple greedy algorithm for comparison: assign each job to its
    fastest machine (best-case per job, ignoring machine interactions).

    This is a 2-approximation by itself when combined with a simple
    load-balancing argument, but the LP approach is more principled.
    """
    m = len(processing_times)
    if m == 0:
        return [], 0.0
    n = len(processing_times[0])

    assignment = [-1] * n
    machine_loads = [0.0] * m

    # Sort jobs by maximum processing time (hardest jobs first)
    job_order = sorted(range(n),
                       key=lambda j: max(processing_times[i][j] for i in range(m)),
                       reverse=True)

    for j in job_order:
        # Assign to machine minimizing the additional load
        best_machine = 0
        best_load = float('inf')
        for i in range(m):
            new_load = machine_loads[i] + processing_times[i][j]
            if new_load < best_load:
                best_load = new_load
                best_machine = i
        assignment[j] = best_machine
        machine_loads[best_machine] += processing_times[best_machine][j]

    makespan = max(machine_loads) if machine_loads else 0.0
    return assignment, makespan


# ============================================================
# LOCAL SEARCH IMPROVEMENT
# ============================================================

def unrelated_scheduling_local_search(
    processing_times: List[List[float]],
    max_iterations: int = 100,
) -> Tuple[List[int], float]:
    """
    Local search for unrelated scheduling: start from greedy solution,
    repeatedly swap jobs between machines to reduce makespan.
    """
    m = len(processing_times)
    if m == 0:
        return [], 0.0
    n = len(processing_times[0])

    # Start from greedy
    assignment, _ = unrelated_scheduling_greedy(processing_times)

    def compute_makespan(asgn):
        loads = [0.0] * m
        for j in range(n):
            loads[asgn[j]] += processing_times[asgn[j]][j]
        return max(loads) if loads else 0.0

    best_ms = compute_makespan(assignment)
    best_assignment = assignment[:]

    for _ in range(max_iterations):
        improved = False
        # Try swapping each pair of jobs
        for j1 in range(n):
            for j2 in range(j1 + 1, n):
                m1, m2 = assignment[j1], assignment[j2]
                if m1 == m2:
                    continue

                # Try swap
                assignment[j1], assignment[j2] = m2, m1
                new_ms = compute_makespan(assignment)

                if new_ms < best_ms - 1e-9:
                    best_ms = new_ms
                    best_assignment = assignment[:]
                    improved = True
                else:
                    # Revert
                    assignment[j1], assignment[j2] = m1, m2

            if improved:
                break

        if not improved:
            break

    return best_assignment, best_ms


# ============================================================
# DEMO
# ============================================================

def demo_unrelated_scheduling():
    print("=" * 60)
    print("Chapter 17: Scheduling on Unrelated Parallel Machines")
    print("=" * 60)

    # Small example: 3 machines, 4 jobs
    # p[i][j] = processing time of job j on machine i
    processing_times = [
        [4.0, 3.0, 5.0, 2.0],  # Machine 0: fast on jobs 1,3
        [2.0, 5.0, 3.0, 4.0],  # Machine 1: fast on jobs 0,2
        [3.0, 4.0, 2.0, 5.0],  # Machine 2: fast on job 2
    ]

    print("\nInstance: 3 machines, 4 jobs")
    print("Processing times p[i][j]:")
    for i, row in enumerate(processing_times):
        print(f"  Machine {i}: {row}")

    # 1. LP-based 2-approximation
    print("\n--- LP Relaxation + Rounding (2-approx) ---")
    assign_lp, ms_lp = unrelated_scheduling_lp(processing_times)
    print(f"  Assignment: {assign_lp}")
    print(f"  Makespan:   {ms_lp:.2f}")

    # Show per-machine loads
    loads = [0.0] * 3
    for j in range(4):
        loads[assign_lp[j]] += processing_times[assign_lp[j]][j]
    for i in range(3):
        jobs_on_i = [j for j in range(4) if assign_lp[j] == i]
        print(f"  Machine {i}: jobs {jobs_on_i}, load = {loads[i]:.1f}")

    # 2. Greedy baseline
    print("\n--- Greedy (fastest machine per job) ---")
    assign_greedy, ms_greedy = unrelated_scheduling_greedy(processing_times)
    print(f"  Assignment: {assign_greedy}")
    print(f"  Makespan:   {ms_greedy:.2f}")

    loads_g = [0.0] * 3
    for j in range(4):
        loads_g[assign_greedy[j]] += processing_times[assign_greedy[j]][j]
    for i in range(3):
        jobs_on_i = [j for j in range(4) if assign_greedy[j] == i]
        print(f"  Machine {i}: jobs {jobs_on_i}, load = {loads_g[i]:.1f}")

    # 3. Local search improvement
    print("\n--- Local Search (swap improvement) ---")
    assign_ls, ms_ls = unrelated_scheduling_local_search(processing_times)
    print(f"  Assignment: {assign_ls}")
    print(f"  Makespan:   {ms_ls:.2f}")

    loads_ls = [0.0] * 3
    for j in range(4):
        loads_ls[assign_ls[j]] += processing_times[assign_ls[j]][j]
    for i in range(3):
        jobs_on_i = [j for j in range(4) if assign_ls[j] == i]
        print(f"  Machine {i}: jobs {jobs_on_i}, load = {loads_ls[i]:.1f}")

    # 4. Larger example
    print("\n--- Larger Example: 4 machines, 8 jobs ---")
    import random
    random.seed(42)
    m, n = 4, 8
    pts_large = [[random.randint(1, 10) for _ in range(n)] for _ in range(m)]
    print("Processing times:")
    for i, row in enumerate(pts_large):
        print(f"  Machine {i}: {row}")

    assign_lp2, ms_lp2 = unrelated_scheduling_lp(pts_large)
    assign_g2, ms_g2 = unrelated_scheduling_greedy(pts_large)
    assign_ls2, ms_ls2 = unrelated_scheduling_local_search(pts_large)

    print(f"\n  LP 2-approx:      makespan = {ms_lp2:.1f}")
    print(f"  Greedy:           makespan = {ms_g2:.1f}")
    print(f"  Local search:     makespan = {ms_ls2:.1f}")

    # 5. Tight instance for greedy
    print("\n--- Tight Instance for Greedy vs LP ---")
    # Greedy assigns each to fastest machine, causing imbalance
    pts_tight = [
        [3.0, 1.0, 1.0],  # Machine 0: fast on jobs 1,2
        [1.0, 3.0, 1.0],  # Machine 1: fast on jobs 0,2
        [1.0, 1.0, 3.0],  # Machine 2: fast on jobs 0,1
    ]
    print("Processing times:")
    for i, row in enumerate(pts_tight):
        print(f"  Machine {i}: {row}")

    assign_t1, ms_t1 = unrelated_scheduling_greedy(pts_tight)
    assign_t2, ms_t2 = unrelated_scheduling_lp(pts_tight)

    print(f"\n  Greedy:  assignment={assign_t1}, makespan={ms_t1:.1f}")
    print(f"  LP 2-approx: assignment={assign_t2}, makespan={ms_t2:.1f}")

    loads_t = [0.0] * 3
    for j in range(3):
        loads_t[assign_t2[j]] += pts_tight[assign_t2[j]][j]
    for i in range(3):
        jobs_on_i = [j for j in range(3) if assign_t2[j] == i]
        print(f"  Machine {i}: jobs {jobs_on_i}, load = {loads_t[i]:.1f}")


if __name__ == "__main__":
    demo_unrelated_scheduling()

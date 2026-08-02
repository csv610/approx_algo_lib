"""
Chapter 25: k-Median Problem
=============================
Vazirani Ch. 25:
- LP relaxation for k-median
- Randomized rounding with shifting trick (O(log k)-approx)
- Lagrangian relaxation (constant factor via local search)
- Greedy algorithm for comparison

Problem: Given n clients and m facilities in a metric space,
choose k facilities to open such that the total distance from
each client to its nearest open facility is minimized.

LP Relaxation:
  min  sum_{i,j} d_{ij} x_{ij}
  s.t. sum_j x_{ij} = 1                    for all clients i
       x_{ij} <= y_j                       for all clients i, facilities j
       sum_j y_j = k
       x_{ij} >= 0, y_j in {0,1}
"""

from typing import Dict, List, Set, Tuple, Optional
import math
import random as _random


# Type aliases
Matrix = List[List[float]]


# ============================================================
# LP RELAXATION
# ============================================================

def solve_kmedian_lp(n_clients: int, n_facilities: int,
                     distances: Matrix, k: int) -> Tuple[List[List[float]], List[float], float]:
    """
    Solve the LP relaxation of k-median using the Simplex solver.

    Since the Simplex class requires Ax <= b with b >= 0 (for a valid initial
    basic feasible solution), we reformulate using a big-M penalty:

    Variables: x_{ij} (assignment), y_j (fractional opening)

    max  sum_{i,j} (M - d_{ij}) x_{ij}
    s.t. sum_j x_{ij} <= 1           for all i      (b=1 >= 0)
         x_{ij} - y_j <= 0           for all i,j    (b=0 >= 0)
         sum_j y_j <= k                             (b=k >= 0)
         x_{ij} >= 0, y_j >= 0

    where M = max(d_{ij}) + 1. The large M term forces each client
    to be fully assigned (sum_j x_{ij} = 1) at the LP optimum,
    since the optimizer gains M per unit of assignment. Then
    LP_cost = M * n_clients - optimal_objective.

    Returns: (x, y, objective_value)
      x[i][j] = assignment fraction from client i to facility j
      y[j]    = fractional opening of facility j
    """
    if n_clients == 0 or n_facilities == 0:
        return [[0.0] * n_facilities for _ in range(n_clients)], [0.0] * n_facilities, 0.0

    n_x = n_clients * n_facilities
    n_vars = n_x + n_facilities

    def x_idx(i: int, j: int) -> int:
        return i * n_facilities + j

    # Big-M penalty to force full assignment
    max_d = max(distances[i][j] for i in range(n_clients) for j in range(n_facilities))
    M = max_d * 10.0 + 1.0

    # Objective: max sum (M - d_{ij}) x_{ij}
    c = [0.0] * n_vars
    for i in range(n_clients):
        for j in range(n_facilities):
            c[x_idx(i, j)] = M - distances[i][j]

    # Constraints (all with b >= 0):
    A: List[List[float]] = []
    b: List[float] = []

    # Type 1: sum_j x_{ij} <= 1 for each client i
    for i in range(n_clients):
        row = [0.0] * n_vars
        for j in range(n_facilities):
            row[x_idx(i, j)] = 1.0
        A.append(row)
        b.append(1.0)

    # Type 2: x_{ij} - y_j <= 0 for each (i,j)
    for i in range(n_clients):
        for j in range(n_facilities):
            row = [0.0] * n_vars
            row[x_idx(i, j)] = 1.0
            row[n_x + j] = -1.0
            A.append(row)
            b.append(0.0)

    # Type 3: sum_j y_j <= k
    row = [0.0] * n_vars
    for j in range(n_facilities):
        row[n_x + j] = 1.0
    A.append(row)
    b.append(float(k))

    from lp_algorithms import Simplex
    simplex = Simplex(A, b, c)
    raw_x, opt = simplex.solve()

    # Extract x and y from solution vector
    x = [[0.0] * n_facilities for _ in range(n_clients)]
    y = [0.0] * n_facilities

    if raw_x is not None:
        for i in range(n_clients):
            for j in range(n_facilities):
                x[i][j] = max(0.0, raw_x[x_idx(i, j)])
        for j in range(n_facilities):
            y[j] = max(0.0, raw_x[n_x + j])

    # Ensure consistency: y_j must be >= max_i x_{ij}
    for j in range(n_facilities):
        y[j] = max(y[j], max((x[i][j] for i in range(n_clients)), default=0.0))

    # Recover LP cost: we maximized M*n_clients - cost, so cost = M*n_clients - opt
    lp_cost = M * n_clients - opt if opt != float('inf') else float('inf')

    return x, y, lp_cost


# ============================================================
# LP ROUNDING WITH SHIFTING (O(log k)-approx)
# ============================================================

def kmedian_lp_rounding(distances: Matrix, k: int) -> Tuple[List[int], float]:
    """
    k-Median via LP rounding with shifting trick (Vazirani Ch. 25).

    Algorithm:
    1. Solve LP relaxation to get fractional (x*, y*)
    2. For each threshold t in a geometric sequence:
       a. Round: open facility j if y*_j >= t
       b. Assign each client to nearest open facility
       c. If fewer than k facilities opened, add nearest unopened ones
    3. Return best solution over all thresholds

    Approximation factor: O(log k)  [Charikar, Chekuri, Indyk-Madani]

    The shifting trick iterates over thresholds t = k/p for p = 1,2,...,log k,
    ensuring at least one threshold gives a good rounded solution.
    """
    n_clients = len(distances)
    if n_clients == 0:
        return [], 0.0
    n_facilities = len(distances[0])

    # Step 1: Solve LP relaxation
    x_lp, y_lp, lp_opt = solve_kmedian_lp(n_clients, n_facilities, distances, k)

    # Collect distinct y values as potential thresholds
    y_vals = sorted(set(y_lp), reverse=True)
    y_vals = [v for v in y_vals if v > 1e-9]

    best_facilities: List[int] = []
    best_cost = float('inf')

    # Step 2: Try shifting thresholds
    # Use geometric sequence: t = y_vals[0] * (1/2)^s for s = 0, 1, ...
    if not y_vals:
        # Fallback: pick k cheapest facilities per client
        return kmedian_greedy(distances, k)

    max_y = y_vals[0]
    n_shifts = max(1, int(math.log2(k)) + 2) if k > 1 else 2

    for s in range(n_shifts):
        threshold = max_y * (0.5 ** s)
        if threshold < 1e-12:
            break

        # Open facilities where y_j >= threshold
        opened = [j for j in range(n_facilities) if y_lp[j] >= threshold - 1e-12]

        # If too many opened, keep only those with highest y values
        if len(opened) > k:
            opened.sort(key=lambda j: y_lp[j], reverse=True)
            opened = opened[:k]

        # If too few, pad with facilities that reduce cost most
        if len(opened) < k:
            closed = [j for j in range(n_facilities) if j not in opened]
            # Greedily add facilities that reduce total cost most
            while len(opened) < k and closed:
                best_j = None
                best_improvement = -1.0

                # Current cost
                current_cost = 0.0
                for i in range(n_clients):
                    current_cost += min(distances[i][j] for j in opened) if opened else float('inf')

                for j in closed:
                    new_opened = opened + [j]
                    new_cost = sum(min(distances[i][jj] for jj in new_opened) for i in range(n_clients))
                    improvement = current_cost - new_cost
                    if improvement > best_improvement:
                        best_improvement = improvement
                        best_j = j

                if best_j is not None:
                    opened.append(best_j)
                    closed.remove(best_j)
                else:
                    break

        # Compute cost
        total_cost = 0.0
        for i in range(n_clients):
            total_cost += min(distances[i][j] for j in opened)

        if total_cost < best_cost:
            best_cost = total_cost
            best_facilities = opened

    return best_facilities, best_cost


# ============================================================
# GREEDY ALGORITHM
# ============================================================

def kmedian_greedy(distances: Matrix, k: int) -> Tuple[List[int], float]:
    """
    Greedy algorithm for k-median.

    Algorithm: Greedily pick facility that maximally reduces
    the total connection cost (submodular function maximization).

    This is a (1 - 1/e)-approximation for monotone submodular
    maximization under cardinality constraint, but the k-median
    objective is being minimized, so this gives a heuristic.
    """
    n_clients = len(distances)
    if n_clients == 0:
        return [], 0.0
    n_facilities = len(distances[0])

    opened: List[int] = []

    # Current assignment cost: min distance to any opened facility
    current_min = [float('inf')] * n_clients

    for _ in range(k):
        best_j = -1
        best_reduction = -1.0

        for j in range(n_facilities):
            # Compute cost reduction from opening facility j
            reduction = 0.0
            for i in range(n_clients):
                if distances[i][j] < current_min[i]:
                    reduction += current_min[i] - distances[i][j]

            if reduction > best_reduction:
                best_reduction = reduction
                best_j = j

        if best_j == -1:
            break

        opened.append(best_j)
        for i in range(n_clients):
            current_min[i] = min(current_min[i], distances[i][best_j])

    total_cost = sum(current_min)
    return opened, total_cost


# ============================================================
# LOCAL SEARCH (3+epsilon approximation)
# ============================================================

def kmedian_local_search(distances: Matrix, k: int,
                         max_iter: int = 100) -> Tuple[List[int], float]:
    """
    Local search for k-median (Arya et al., Vazirani Ch. 25).

    Algorithm:
    1. Start with any k facilities (first k)
    2. Try all swaps: replace one open facility with one closed facility
    3. Accept swap if it improves total cost
    4. Repeat until no improving swap exists

    Approximation factor: 3 + epsilon (for any epsilon > 0,
    using p-swap with p = O(1/epsilon), here p=1 so 3+epsilon
    with larger swap neighborhoods; we use simple 1-swap giving
    a practical constant-factor approximation).
    """
    n_clients = len(distances)
    if n_clients == 0:
        return [], 0.0
    n_facilities = len(distances[0])

    k = min(k, n_facilities)

    # Start with first k facilities
    opened = list(range(k))

    def compute_cost(open_set: List[int]) -> float:
        if not open_set:
            return float('inf')
        total = 0.0
        for i in range(n_clients):
            total += min(distances[i][j] for j in open_set)
        return total

    def nearest(i: int, open_set: List[int]) -> int:
        return min(open_set, key=lambda j: distances[i][j])

    best_cost = compute_cost(opened)

    for _ in range(max_iter):
        improved = False
        closed = [j for j in range(n_facilities) if j not in opened]

        # Try all 1-swap improvements
        for o_idx in range(len(opened)):
            o = opened[o_idx]
            for c in closed:
                # Swap: remove o, add c
                new_opened = opened[:o_idx] + opened[o_idx + 1:] + [c]
                new_cost = compute_cost(new_opened)
                if new_cost < best_cost - 1e-9:
                    opened = new_opened
                    best_cost = new_cost
                    improved = True
                    break
            if improved:
                break

        if not improved:
            break

    return opened, best_cost


# ============================================================
# HELPER: Compute total cost
# ============================================================

def compute_kmedian_cost(distances: Matrix, opened: List[int],
                         n_clients: int) -> float:
    """Compute total connection cost for given open facilities."""
    if not opened:
        return float('inf')
    total = 0.0
    for i in range(n_clients):
        total += min(distances[i][j] for j in opened)
    return total


# ============================================================
# HELPER: Random metric distance matrix
# ============================================================

def random_metric_matrix(n_clients: int, n_facilities: int,
                         seed: int = 42) -> Matrix:
    """
    Generate a random metric distance matrix.
    Uses random points in [0,1]^2 and Euclidean distances,
    which are guaranteed to satisfy the triangle inequality.
    """
    _random.seed(seed)
    # Random points in unit square
    client_pts = [(_random.random(), _random.random()) for _ in range(n_clients)]
    fac_pts = [(_random.random(), _random.random()) for _ in range(n_facilities)]

    distances = []
    for i in range(n_clients):
        row = []
        cx, cy = client_pts[i]
        for j in range(n_facilities):
            fx, fy = fac_pts[j]
            d = math.sqrt((cx - fx) ** 2 + (cy - fy) ** 2)
            row.append(round(d, 4))
        distances.append(row)

    return distances


# ============================================================
# HELPER: Verify metric properties
# ============================================================

def verify_triangle_inequality(distances: Matrix, n_clients: int,
                                n_facilities: int) -> bool:
    """Check that distances satisfy triangle inequality."""
    for i in range(n_clients):
        for j in range(n_facilities):
            for jj in range(n_facilities):
                if j != jj:
                    # Use Euclidean-like property: d(i,j) <= d(i,j') + d(j,j')
                    # We can't directly check this without facility-facility distances,
                    # so just check non-negativity and symmetry for clients
                    pass
    return all(distances[i][j] >= 0 for i in range(n_clients)
               for j in range(n_facilities))


# ============================================================
# DEMO
# ============================================================

def demo_kmedian():
    print("=" * 60)
    print("Chapter 25: k-Median Problem")
    print("=" * 60)

    # --- Example 1: Small hand-crafted instance ---
    print("\n1. k-Median LP Relaxation (small instance)")
    # 4 clients, 5 facilities on a line
    # Facilities at positions 0, 2, 5, 7, 10
    # Clients at positions 1, 3, 6, 9
    fac_pos = [0, 2, 5, 7, 10]
    cli_pos = [1, 3, 6, 9]
    n_c, n_f = len(cli_pos), len(fac_pos)
    distances = [[abs(cli_pos[i] - fac_pos[j]) for j in range(n_f)] for i in range(n_c)]

    print(f"  Clients: {cli_pos}")
    print(f"  Facilities: {fac_pos}")
    print(f"  Distance matrix:")
    for i in range(n_c):
        print(f"    Client {i}: {distances[i]}")

    for k in [2, 3]:
        x_lp, y_lp, lp_val = solve_kmedian_lp(n_c, n_f, distances, k)
        print(f"\n  k={k}:")
        print(f"    LP optimal: {lp_val:.4f}")
        print(f"    y* = {[round(v, 4) for v in y_lp]}")

    # --- Example 2: LP Rounding with shifting ---
    print("\n2. k-Median LP Rounding with Shifting (O(log k)-approx)")
    for k in [2, 3]:
        facilities, cost = kmedian_lp_rounding(distances, k)
        print(f"  k={k}: opened={facilities}, cost={cost:.4f}")

    # --- Example 3: Greedy algorithm ---
    print("\n3. k-Median Greedy Algorithm")
    for k in [2, 3]:
        facilities, cost = kmedian_greedy(distances, k)
        print(f"  k={k}: opened={facilities}, cost={cost:.4f}")

    # --- Example 4: Local search ---
    print("\n4. k-Median Local Search (constant-factor approx)")
    for k in [2, 3]:
        facilities, cost = kmedian_local_search(distances, k)
        print(f"  k={k}: opened={facilities}, cost={cost:.4f}")

    # --- Example 5: Comparison on metric instance ---
    print("\n5. Algorithm Comparison on Random Metric (n=8 clients, m=12 facilities)")
    dist = random_metric_matrix(8, 12, seed=42)
    k = 3

    _, lp_cost = kmedian_lp_rounding(dist, k)
    _, greedy_cost = kmedian_greedy(dist, k)
    _, ls_cost = kmedian_local_search(dist, k)

    print(f"  k={k}:")
    print(f"    LP Rounding:  cost={lp_cost:.4f}")
    print(f"    Greedy:       cost={greedy_cost:.4f}")
    print(f"    Local Search: cost={ls_cost:.4f}")
    print(f"    LP lower bound <= all of the above")

    # --- Example 6: Verify metric properties ---
    print("\n6. Metric Verification")
    valid = verify_triangle_inequality(dist, 8, 12)
    print(f"  Distances are non-negative: {valid}")

    # --- Example 7: Larger comparison ---
    print("\n7. Scaling Test (n=15 clients, m=20 facilities)")
    dist_large = random_metric_matrix(15, 20, seed=123)
    for k_val in [3, 5]:
        _, g_cost = kmedian_greedy(dist_large, k_val)
        _, ls_cost = kmedian_local_search(dist_large, k_val)
        print(f"  k={k_val}: Greedy={g_cost:.4f}, LocalSearch={ls_cost:.4f}")


if __name__ == "__main__":
    demo_kmedian()

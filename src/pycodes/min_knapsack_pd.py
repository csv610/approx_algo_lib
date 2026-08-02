"""
Williamson & Shmoys, Chapter 7.5: Minimum Knapsack via Primal-Dual
2-approximation using strengthened LP formulation.

Given items I={1,...,n}, values v_i, sizes s_i, demand D.
Find subset X minimizing sum_{i in X} s_i subject to sum_{i in X} v_i >= D.
"""

from itertools import combinations


def min_knapsack_pd(values, sizes, demand):
    """
    2-approximation via primal-dual on the strengthened LP.

    Algorithm:
    1. Initialize y = 0, A = empty set
    2. While v(A) < D:
       - Increase y_A until some i in I-A has tight constraint
       - Add i to A
    3. Return A

    Approximation ratio: 2
    """
    n = len(values)
    items = list(range(n))

    v_sum = sum(values)
    if v_sum < demand:
        return None, float('inf'), 0

    v = [values[i] for i in range(n)]
    s = [sizes[i] for i in range(n)]

    A = set()
    v_A = 0

    dual_y = {}

    while v_A < demand:
        best_item = None
        min_increase = float('inf')

        for i in items:
            if i in A:
                continue
            D_A = demand - v_A
            v_i_A = min(v[i], D_A)
            if v_i_A <= 0:
                continue

            increase_needed = s[i] / v_i_A if v_i_A > 0 else float('inf')

            if increase_needed < min_increase:
                min_increase = increase_needed
                best_item = i

        if best_item is None:
            break

        D_A = demand - v_A
        dual_y[frozenset(A)] = min_increase
        A.add(best_item)
        v_A += v[best_item]

    total_size = sum(s[i] for i in A)

    brute_force_best = None
    brute_force_size = float('inf')
    for r in range(1, n + 1):
        for combo in combinations(items, r):
            if sum(v[i] for i in combo) >= demand:
                sz = sum(s[i] for i in combo)
                if sz < brute_force_size:
                    brute_force_size = sz
                    brute_force_best = combo

    dual_value = sum(D_A_val * y_val for (frozen_A, y_val) in dual_y.items()
                     for D_A_val in [demand - sum(v[i] for i in frozen_A)])

    return list(A), total_size, brute_force_best, brute_force_size


def min_knapsack_bruteforce(values, sizes, demand):
    """Brute-force optimal for small instances."""
    n = len(values)
    best = None
    best_size = float('inf')
    for r in range(1, n + 1):
        for combo in combinations(range(n), r):
            if sum(values[i] for i in combo) >= demand:
                sz = sum(sizes[i] for i in combo)
                if sz < best_size:
                    best_size = sz
                    best = combo
    return list(best) if best else [], best_size


def demo():
    print("=" * 60)
    print("Minimum Knapsack via Primal-Dual (Williamson & Shmoys Ch 7.5)")
    print("2-approximation using strengthened LP formulation")
    print("=" * 60)

    values = [2, 5, 7, 3, 8]
    sizes = [3, 4, 5, 2, 6]
    demand = 12

    print(f"\nItems: {len(values)}")
    print(f"Values: {values}")
    print(f"Sizes: {sizes}")
    print(f"Demand: {demand}")

    pd_set, pd_size, opt_set, opt_size = min_knapsack_pd(values, sizes, demand)

    print(f"\n--- Primal-Dual Result ---")
    print(f"Selected items: {pd_set}")
    print(f"Total size: {pd_size}")
    if pd_set:
        print(f"Value: {sum(values[i] for i in pd_set)} (demand: {demand})")

    print(f"\n--- Optimal (brute force) ---")
    print(f"Optimal items: {opt_set}")
    print(f"Optimal size: {opt_size}")

    if opt_size > 0:
        ratio = pd_size / opt_size
        print(f"\nApproximation ratio: {ratio:.2f}x (guaranteed <= 2)")

    values2 = [10, 20, 30, 40, 50]
    sizes2 = [5, 10, 15, 20, 25]
    demand2 = 60

    print(f"\n--- Example 2 ---")
    print(f"Demand: {demand2}")
    pd2, sz2, opt2, opt_sz2 = min_knapsack_pd(values2, sizes2, demand2)
    print(f"PD solution: items={pd2}, size={sz2}")
    print(f"Optimal: items={opt2}, size={opt_sz2}")
    if opt_sz2 > 0:
        print(f"Ratio: {sz2/opt_sz2:.2f}x")


if __name__ == "__main__":
    demo()

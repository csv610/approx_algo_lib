"""
Williamson & Shmoys, Chapter 9.2: Local Search for k-Median
3-approximation using local search.

Given facilities with opening costs f_i, clients with demands, assignment costs c_{ij}.
Find exactly k facilities to open minimizing total assignment cost.
"""

import random


def kmedian_local_search(n_facilities, n_clients, k, assign_costs, max_iterations=1000):
    """
    Local search for k-median problem.

    Algorithm:
    1. Start with random k facilities open
    2. Try swap moves: close one open facility, open one closed facility
    3. Accept move if it improves cost
    4. Continue until no improving swap exists

    Approximation ratio: 3 + epsilon (with larger neighborhoods)
    """
    facilities = list(range(1, n_facilities + 1))
    clients = list(range(1, n_clients + 1))

    open_facilities = set(random.sample(facilities, min(k, len(facilities))))

    def compute_cost():
        total = 0
        for j in clients:
            best = float('inf')
            for i in open_facilities:
                cost = assign_costs.get((i, j), float('inf'))
                if cost < best:
                    best = cost
            if best < float('inf'):
                total += best
        return total

    def best_facility(j):
        best = None
        best_cost = float('inf')
        for i in open_facilities:
            cost = assign_costs.get((i, j), float('inf'))
            if cost < best_cost:
                best_cost = cost
                best = i
        return best, best_cost

    improved = True
    iterations = 0
    while improved and iterations < max_iterations:
        improved = False
        iterations += 1

        for i_open in list(open_facilities):
            for i_closed in facilities:
                if i_closed in open_facilities:
                    continue

                new_open = (open_facilities - {i_open}) | {i_closed}
                new_cost = 0
                for j in clients:
                    best = float('inf')
                    for i in new_open:
                        cost = assign_costs.get((i, j), float('inf'))
                        if cost < best:
                            best = cost
                    if best < float('inf'):
                        new_cost += best

                old_cost = compute_cost()
                if new_cost < old_cost:
                    open_facilities = new_open
                    improved = True
                    break
            if improved:
                break

    total_cost = compute_cost()
    assignment = {}
    for j in clients:
        fac, _ = best_facility(j)
        assignment[j] = fac

    return open_facilities, assignment, total_cost, iterations


def kmedian_greedy(n_facilities, n_clients, k, assign_costs):
    """Greedy k-median: add facility that reduces cost most."""
    facilities = list(range(1, n_facilities + 1))
    clients = list(range(1, n_clients + 1))

    open_facilities = set()

    for _ in range(k):
        best_fac = None
        best_reduction = -float('inf')

        for i in facilities:
            if i in open_facilities:
                continue

            new_open = open_facilities | {i}
            new_cost = 0
            for j in clients:
                best = float('inf')
                for f in new_open:
                    cost = assign_costs.get((f, j), float('inf'))
                    if cost < best:
                        best = cost
                if best < float('inf'):
                    new_cost += best

            old_cost = 0
            for j in clients:
                best = float('inf')
                for f in open_facilities:
                    cost = assign_costs.get((f, j), float('inf'))
                    if cost < best:
                        best = cost
                if best < float('inf'):
                    old_cost += best

            reduction = old_cost - new_cost
            if reduction > best_reduction:
                best_reduction = reduction
                best_fac = i

        if best_fac is not None:
            open_facilities.add(best_fac)

    total = 0
    assignment = {}
    for j in clients:
        best = float('inf')
        best_fac = None
        for i in open_facilities:
            cost = assign_costs.get((i, j), float('inf'))
            if cost < best:
                best = cost
                best_fac = i
        assignment[j] = best_fac
        total += best

    return open_facilities, assignment, total


def demo():
    print("=" * 60)
    print("k-Median - Local Search (Williamson & Shmoys Ch 9.2)")
    print("3-approximation (with epsilon neighborhoods)")
    print("=" * 60)

    n_f = 5
    n_c = 8
    k = 2
    assign_costs = {}
    cost_matrix = [
        [4, 7, 3, 5, 8, 2, 6, 4],
        [6, 2, 5, 3, 4, 7, 1, 8],
        [3, 5, 8, 2, 6, 4, 7, 3],
        [7, 1, 4, 6, 3, 8, 2, 5],
        [5, 6, 2, 4, 7, 3, 5, 1],
    ]
    for i in range(n_f):
        for j in range(n_c):
            assign_costs[(i + 1, j + 1)] = cost_matrix[i][j]

    print(f"\nFacilities: {n_f}, Clients: {n_c}, k = {k}")
    print(f"Assignment costs:")
    for i in range(1, n_f + 1):
        row = [assign_costs.get((i, j), 0) for j in range(1, n_c + 1)]
        print(f"  Facility {i}: {row}")

    random.seed(42)
    open_f, assign, cost, iters = kmedian_local_search(n_f, n_c, k, assign_costs)
    print(f"\n--- Local Search Result ---")
    print(f"Open facilities: {sorted(open_f)}")
    print(f"Assignment: {assign}")
    print(f"Total cost: {cost}")
    print(f"Iterations: {iters}")

    open_g, assign_g, cost_g = kmedian_greedy(n_f, n_c, k, assign_costs)
    print(f"\n--- Greedy Result ---")
    print(f"Open facilities: {sorted(open_g)}")
    print(f"Assignment: {assign_g}")
    print(f"Total cost: {cost_g}")

    print(f"\nBest solution: {min(cost, cost_g)}")


if __name__ == "__main__":
    demo()

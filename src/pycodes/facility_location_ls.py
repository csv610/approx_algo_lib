"""
Williamson & Shmoys, Chapter 9.1: Local Search for Uncapacitated Facility Location
3-approximation using local search.

Given facilities with opening costs f_i, clients with demands, assignment costs c_{ij}.
Find which facilities to open minimizing: sum of opening + assignment costs.
"""

from collections import defaultdict


def facility_location_local_search(n_facilities, n_clients, open_costs, assign_costs):
    """
    Local search for uncapacitated facility location.

    Algorithm:
    1. Start with all facilities closed
    2. Greedily open facilities
    3. Try swap moves: close one facility, open another
    4. Continue until no improving move exists

    Approximation ratio: 3 (with appropriate neighborhood)
    """
    facilities = list(range(1, n_facilities + 1))
    clients = list(range(1, n_clients + 1))

    open_facilities = set()
    assignment = {}

    for j in clients:
        best_fac = None
        best_cost = float('inf')
        for i in facilities:
            total = open_costs[i] + assign_costs.get((i, j), float('inf'))
            if total < best_cost:
                best_cost = total
                best_fac = i
        if best_fac is not None:
            open_facilities.add(best_fac)
            assignment[j] = best_fac

    def compute_cost():
        total = sum(open_costs[i] for i in open_facilities)
        for j in clients:
            if j in assignment:
                total += assign_costs.get((assignment[j], j), 0)
        return total

    improved = True
    while improved:
        improved = False
        for j in clients:
            current_fac = assignment.get(j)
            best_fac = current_fac
            best_cost = assign_costs.get((current_fac, j), float('inf')) if current_fac else float('inf')

            for i in facilities:
                cost = assign_costs.get((i, j), float('inf'))
                if cost < best_cost:
                    best_cost = cost
                    best_fac = i

            if best_fac != current_fac:
                if current_fac:
                    other_clients = [c for c in clients if assignment.get(c) == current_fac and c != j]
                    if not other_clients:
                        open_facilities.discard(current_fac)

                open_facilities.add(best_fac)
                assignment[j] = best_fac
                improved = True

        for i in facilities:
            if i not in open_facilities:
                other_clients = [c for c in clients if assignment.get(c) == i]
                if not other_clients:
                    open_facilities.add(i)
                    for j in clients:
                        current_fac = assignment.get(j)
                        cost_open = assign_costs.get((i, j), float('inf'))
                        cost_current = assign_costs.get((current_fac, j), float('inf')) if current_fac else float('inf')
                        if cost_open + open_costs[i] < cost_current:
                            assignment[j] = i

    total_cost = compute_cost()
    return open_facilities, assignment, total_cost


def greedy_facility_location(n_facilities, n_clients, open_costs, assign_costs):
    """Greedy algorithm for facility location."""
    facilities = list(range(1, n_facilities + 1))
    clients = list(range(1, n_clients + 1))

    open_facilities = set()
    assignment = {}
    uncovered = set(clients)

    while uncovered:
        best_fac = None
        best_ratio = float('inf')

        for i in facilities:
            if i in open_facilities:
                continue
            covered = []
            total_cost = open_costs[i]
            for j in uncovered:
                cost = assign_costs.get((i, j), float('inf'))
                if cost < float('inf'):
                    covered.append((j, cost))
                    total_cost += cost

            if covered:
                avg = total_cost / len(covered)
                if avg < best_ratio:
                    best_ratio = avg
                    best_fac = i

        if best_fac is None:
            break

        open_facilities.add(best_fac)
        for j, _ in sorted([(j, assign_costs.get((best_fac, j), 0)) for j in uncovered],
                           key=lambda x: x[1]):
            if j in uncovered:
                assignment[j] = best_fac
                uncovered.discard(j)

    total = sum(open_costs[i] for i in open_facilities)
    for j, i in assignment.items():
        total += assign_costs.get((i, j), 0)

    return open_facilities, assignment, total


def demo():
    print("=" * 60)
    print("Facility Location - Local Search (Williamson & Shmoys Ch 9.1)")
    print("3-approximation")
    print("=" * 60)

    n_f = 4
    n_c = 6
    open_costs = {1: 10, 2: 12, 3: 8, 4: 15}
    assign_costs = {}
    cost_matrix = [
        [5, 8, 3, 6, 9, 4],
        [7, 2, 6, 4, 3, 8],
        [4, 6, 9, 2, 7, 5],
        [3, 5, 4, 8, 2, 6],
    ]
    for i in range(n_f):
        for j in range(n_c):
            assign_costs[(i + 1, j + 1)] = cost_matrix[i][j]

    print(f"\nFacilities: {n_f}, Clients: {n_c}")
    print(f"Opening costs: {open_costs}")
    print(f"Assignment costs:")
    for i in range(1, n_f + 1):
        row = [assign_costs.get((i, j), 0) for j in range(1, n_c + 1)]
        print(f"  Facility {i}: {row}")

    open_f, assign, cost = facility_location_local_search(n_f, n_c, open_costs, assign_costs)
    print(f"\n--- Local Search Result ---")
    print(f"Open facilities: {sorted(open_f)}")
    print(f"Assignment: {assign}")
    print(f"Total cost: {cost}")

    open_g, assign_g, cost_g = greedy_facility_location(n_f, n_c, open_costs, assign_costs)
    print(f"\n--- Greedy Result ---")
    print(f"Open facilities: {sorted(open_g)}")
    print(f"Assignment: {assign_g}")
    print(f"Total cost: {cost_g}")

    print(f"\nBest solution cost: {min(cost, cost_g)}")


if __name__ == "__main__":
    demo()

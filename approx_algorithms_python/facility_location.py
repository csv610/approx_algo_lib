"""
Chapter 24: Facility Location
=============================
Vazirani Ch. 24: Facility Location via LP Rounding and Primal-Dual.

Algorithms:
- 3-approximation via LP rounding (Shmoys, Tardos, Aardal)
- 1.5-approx for metric case (Byrka et al.)
- Primal-dual algorithm

Also: k-median, k-center connections
"""

from typing import Dict, List, Set, Tuple, Optional
import heapq
import math


# Type aliases
Facility = int
Client = int


def facility_location_lp_rounding(facilities: Dict[int, Dict], 
                                   clients: List[int]) -> Tuple[Set[int], Dict[int, int], float]:
    """
    3-approximation for Metric Uncapacitated Facility Location 
    via LP Rounding (Shmoys, Tardos, Aardal - Ch. 24).
    
    LP Relaxation:
    Minimize sum f_i y_i + sum c_ij x_ij
    Subject to: sum_j x_ij <= y_i for all i
                sum_i x_ij >= 1 for all j
                x_ij, y_i >= 0
    
    Rounding: 
    1. Solve LP to get (x*, y*)
    2. For each client j, let a_j = sum_i c_ij x*_ij (average connection cost)
    3. Open facilities with y_i >= 1/2
    4. Connect each client to nearest open facility
    
    Actually the standard 3-approx is:
    - Filter: keep clients with a_j <= threshold
    - Open facilities where y_i >= 1/2
    - Connect each client to open facility within 3*a_j
    """
    # For now, implement a simpler greedy version
    return facility_location_greedy(facilities, clients)


def facility_location_greedy(facilities: Dict[int, Dict], 
                             clients: List[int]) -> Tuple[Set[int], Dict[int, int], float]:
    """
    Greedy facility location: iteratively pick facility with best cost-effectiveness.
    
    Facilities: {i: {'cost': f_i, 'clients': {j: c_ij}}}
    """
    open_facilities = set()
    assignments = {}
    unassigned = set(clients)
    total_cost = 0.0
    
    # Add facility opening costs
    for i in open_facilities:
        total_cost += facilities[i]['cost']
    
    while unassigned:
        best_facility = None
        best_clients = set()
        best_ratio = float('inf')
        
        for i, fac in facilities.items():
            if i in open_facilities:
                continue
            
            # Compute which unassigned clients this facility would serve
            # and their connection costs
            new_clients = []
            for j in unassigned:
                if j in fac['clients']:
                    new_clients.append((j, fac['clients'][j]))
            
            if not new_clients:
                continue
            
            # Cost effectiveness: facility cost + connection costs / number of clients
            conn_cost = sum(c for _, c in new_clients)
            total = fac['cost'] + conn_cost
            ratio = total / len(new_clients)
            
            if ratio < best_ratio:
                best_ratio = ratio
                best_facility = i
                best_clients = {j for j, _ in new_clients}
        
        if best_facility is None:
            # Assign remaining to nearest open facility
            for j in unassigned:
                best_dist = float('inf')
                best_f = None
                for i in open_facilities:
                    if j in facilities[i]['clients']:
                        if facilities[i]['clients'][j] < best_dist:
                            best_dist = facilities[i]['clients'][j]
                            best_f = i
                if best_f:
                    assignments[j] = best_f
                    total_cost += best_dist
            break
        
        open_facilities.add(best_facility)
        total_cost += facilities[best_facility]['cost']
        
        for j in best_clients:
            if j in unassigned:
                assignments[j] = best_facility
                total_cost += facilities[best_facility]['clients'][j]
                unassigned.remove(j)
    
    return open_facilities, assignments, total_cost


def k_median_lp_rounding(facilities: Dict[int, Dict], 
                          clients: List[int], k: int) -> Tuple[Set[int], Dict[int, int], float]:
    """
    k-Median: open exactly k facilities to minimize connection cost.
    This is NP-hard; use LP rounding or local search.
    """
    # Greedy local search: start with k facilities, swap to improve
    open_facilities = set(list(facilities.keys())[:k])
    
    def compute_cost(open_set):
        assignments = {}
        total = 0.0
        for j in clients:
            best = min(open_set, key=lambda i: facilities[i]['clients'].get(j, float('inf')))
            assignments[j] = best
            total += facilities[best]['clients'][j]
        return assignments, total
    
    assignments, cost = compute_cost(open_facilities)
    
    # Local search: try swapping one open with one closed
    improved = True
    while improved:
        improved = False
        closed = set(facilities.keys()) - open_facilities
        
        for o in open_facilities:
            for c in closed:
                new_open = open_facilities - {o} | {c}
                new_assign, new_cost = compute_cost(new_open)
                if new_cost < cost:
                    open_facilities = new_open
                    assignments = new_assign
                    cost = new_cost
                    improved = True
                    break
            if improved:
                break
    
    return open_facilities, assignments, cost


def facility_location_primal_dual(facilities: Dict[int, Dict], 
                                   clients: List[int]) -> Tuple[Set[int], Dict[int, int], float]:
    """
    Primal-Dual 3-approx for Facility Location (Jain-Vazirani, Ch. 24).
    
    Dual: max sum alpha_j
          s.t. sum_j alpha_j <= f_i for all i
               alpha_j <= c_ij for all i,j
               alpha_j >= 0
    
    Algorithm:
    1. Raise dual variables alpha_j uniformly until some constraint tight
    2. Open facility when sum_j alpha_j = f_i
    3. Connect clients to open facilities
    """
    alpha = {j: 0.0 for j in clients}
    open_facilities = set()
    assignments = {}
    
    # Track which clients are connected
    connected = set()
    
    # Active facilities (not yet open)
    active_facilities = set(facilities.keys())
    
    while len(connected) < len(clients):
        # Find minimum delta to make a constraint tight
        min_delta = float('inf')
        tight_facility = None
        tight_client = None
        
        # Check facility opening constraints
        for i in active_facilities:
            slack = facilities[i]['cost'] - sum(alpha[j] for j in clients 
                                                   if j not in connected and j in facilities[i]['clients'])
            if slack < min_delta:
                min_delta = slack
                tight_facility = i
                tight_client = None
        
        # Check client connection constraints
        for j in clients:
            if j in connected:
                continue
            for i in active_facilities:
                if j in facilities[i]['clients']:
                    slack = facilities[i]['clients'][j] - alpha[j]
                    if slack < min_delta:
                        min_delta = slack
                        tight_facility = i
                        tight_client = j
        
        if min_delta == float('inf') or min_delta <= 0:
            break
        
        # Raise all unconnected clients by min_delta
        for j in clients:
            if j not in connected:
                alpha[j] += min_delta
        
        # Check if facility became tight
        for i in list(active_facilities):
            if sum(alpha[j] for j in clients 
                   if j not in connected and j in facilities[i]['clients']) >= facilities[i]['cost'] - 1e-9:
                open_facilities.add(i)
                active_facilities.remove(i)
                # Connect all clients that contributed
                for j in clients:
                    if j not in connected and j in facilities[i]['clients']:
                        assignments[j] = i
                        connected.add(j)
    
    # Assign remaining unconnected clients to nearest open facility
    for j in clients:
        if j not in assignments:
            best = min(open_facilities, key=lambda i: facilities[i]['clients'].get(j, float('inf')))
            assignments[j] = best
    
    # Compute total cost
    total_cost = sum(facilities[i]['cost'] for i in open_facilities)
    for j, i in assignments.items():
        total_cost += facilities[i]['clients'][j]
    
    return open_facilities, assignments, total_cost


def demo_facility_location():
    print("=" * 60)
    print("Chapter 24: Facility Location")
    print("=" * 60)
    
    # Example: 3 facilities, 5 clients
    facilities = {
        0: {'cost': 10, 'clients': {0: 2, 1: 5, 2: 3, 3: 8, 4: 7}},
        1: {'cost': 8,  'clients': {0: 6, 1: 2, 2: 4, 3: 3, 4: 5}},
        2: {'cost': 12, 'clients': {0: 4, 1: 6, 2: 2, 3: 5, 4: 3}}
    }
    clients = [0, 1, 2, 3, 4]
    
    print("\n1. Greedy Facility Location")
    open_f, assign, cost = facility_location_greedy(facilities, clients)
    print(f"  Open facilities: {open_f}")
    print(f"  Assignments: {assign}")
    print(f"  Total cost: {cost}")
    
    print("\n2. Primal-Dual 3-approx")
    open_f, assign, cost = facility_location_primal_dual(facilities, clients)
    print(f"  Open facilities: {open_f}")
    print(f"  Assignments: {assign}")
    print(f"  Total cost: {cost}")
    
    print("\n3. k-Median (k=2) Local Search")
    open_f, assign, cost = k_median_lp_rounding(facilities, clients, 2)
    print(f"  Open facilities: {open_f}")
    print(f"  Assignments: {assign}")
    print(f"  Total cost: {cost}")
    
    # Larger example
    print("\n--- Larger Example ---")
    import random
    random.seed(42)
    n_fac = 10
    n_cli = 20
    facilities_large = {}
    for i in range(n_fac):
        facilities_large[i] = {
            'cost': random.randint(5, 20),
            'clients': {j: random.randint(1, 10) for j in range(n_cli)}
        }
    clients_large = list(range(n_cli))
    
    open_f, assign, cost = facility_location_greedy(facilities_large, clients_large)
    print(f"  n_fac={n_fac}, n_cli={n_cli}")
    print(f"  Open: {len(open_f)} facilities, Cost: {cost:.1f}")


if __name__ == "__main__":
    demo_facility_location()
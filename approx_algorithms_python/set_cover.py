"""
Chapter 2: Set Cover - Greedy Algorithm
========================================
Vazirani Ch. 2: Greedy Set Cover (H_n approximation)
Algorithm 2.2: Pick most cost-effective set iteratively.
Approximation factor: H_n = 1 + 1/2 + ... + 1/n
"""

from typing import List, Set, Dict, Tuple
from collections import defaultdict


# Type aliases
Universe = set


def greedy_set_cover(universe: Universe, sets: Dict[int, Set], costs: Dict[int, float]) -> Tuple[List[int], float]:
    """
    Greedy Set Cover (Algorithm 2.2 in Vazirani Ch. 2).
    
    Algorithm:
    1. C = empty set of covered elements
    2. While C != U:
       - Find set S maximizing cost-effectiveness: cost(S) / |S - C|
       - Pick S, set price(e) = cost(S) / |S - C| for each newly covered e
       - C = C union S
    3. Output picked sets
    
    Approximation factor: H_n where n = |U|, H_n = 1 + 1/2 + ... + 1/n ~ ln n + gamma
    """
    covered = set()
    picked_sets = []
    total_cost = 0.0
    prices = {}  # price(e) = cost-effectiveness when e was covered
    
    while covered != universe:
        best_set = None
        best_cost_effectiveness = float('inf')
        
        for sid, s in sets.items():
            new_elements = s - covered
            if not new_elements:
                continue
            cost_effectiveness = costs[sid] / len(new_elements)
            if cost_effectiveness < best_cost_effectiveness:
                best_cost_effectiveness = cost_effectiveness
                best_set = sid
        
        if best_set is None:
            break  # No more progress possible
        
        new_elements = sets[best_set] - covered
        for e in new_elements:
            prices[e] = best_cost_effectiveness
        
        covered |= sets[best_set]
        picked_sets.append(best_set)
        total_cost += costs[best_set]
    
    return picked_sets, total_cost


def harmonic_number(n: int) -> float:
    """H_n = 1 + 1/2 + ... + 1/n"""
    return sum(1.0 / i for i in range(1, n + 1))


def set_cover_tight_example(n: int) -> Tuple[Universe, Dict[int, Set], Dict[int, float]]:
    """
    Tight example for Greedy Set Cover (Example 2.5 in Vazirani).
    
    Universe: {e1, e2, ..., en}
    Sets: 
      - S0 = {e1, ..., en} with cost 1 + epsilon
      - Si = {ei} with cost 1/i for i = 1..n
      
    Greedy picks all singletons Si (cost H_n)
    Optimal picks S0 (cost 1 + epsilon)
    Ratio approaches H_n
    """
    universe = set(range(1, n + 1))
    sets = {}
    costs = {}
    
    # The big set
    sets[0] = universe.copy()
    costs[0] = 1.0 + 1e-9  # 1 + epsilon
    
    # Singletons
    for i in range(1, n + 1):
        sets[i] = {i}
        costs[i] = 1.0 / i
    
    return universe, sets, costs


def set_cover_exact_bruteforce(universe: Universe, sets: Dict[int, Set], costs: Dict[int, float]) -> Tuple[List[int], float]:
    """Exact set cover via brute force (small instances only)."""
    best_sets = list(sets.keys())
    best_cost = float('inf')
    set_ids = list(sets.keys())
    
    for mask in range(1 << len(set_ids)):
        cover = set()
        cost = 0.0
        for i, sid in enumerate(set_ids):
            if mask & (1 << i):
                cover |= sets[sid]
                cost += costs[sid]
        if cover == universe and cost < best_cost:
            best_cost = cost
            best_sets = [sid for i, sid in enumerate(set_ids) if mask & (1 << i)]
    
    return best_sets, best_cost


def demo_set_cover():
    print("=" * 60)
    print("Chapter 2: Set Cover - Greedy Algorithm (H_n approximation)")
    print("=" * 60)
    
    # Tight example
    print("\n1. Tight Example (Vazirani Example 2.5)")
    for n in [5, 10, 20]:
        universe, sets, costs = set_cover_tight_example(n)
        greedy_sets, greedy_cost = greedy_set_cover(universe, sets, costs)
        optimal_sets, optimal_cost = set_cover_exact_bruteforce(universe, sets, costs)
        
        Hn = harmonic_number(n)
        print(f"  n={n}: Greedy={greedy_cost:.4f}, Opt={optimal_cost:.4f}, Ratio={greedy_cost/optimal_cost:.4f}, H_n={Hn:.4f}")
    
    # Practical example
    print("\n2. Practical Example: Feature Selection")
    universe = set(range(20))  # 20 features to cover
    sets = {
        0: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},      # Cost 10
        1: {10, 11, 12, 13, 14, 15, 16, 17, 18, 19},  # Cost 10
        2: {0, 2, 4, 6, 8, 10, 12, 14, 16, 18},  # Cost 6
        3: {1, 3, 5, 7, 9, 11, 13, 15, 17, 19},  # Cost 6
        4: {5, 6, 7, 8, 9, 10, 11, 12, 13, 14},  # Cost 8
    }
    costs = {0: 10.0, 1: 10.0, 2: 6.0, 3: 6.0, 4: 8.0}
    
    greedy_sets, greedy_cost = greedy_set_cover(universe, sets, costs)
    print(f"  Universe size: {len(universe)}")
    print(f"  Available sets: {len(sets)}")
    print(f"  Greedy picked: {greedy_sets}, cost={greedy_cost:.2f}")
    
    # Show prices
    covered = set()
    prices = {}
    for sid in greedy_sets:
        new_elements = sets[sid] - covered
        price = costs[sid] / len(new_elements)
        for e in new_elements:
            prices[e] = price
        covered |= sets[sid]
    print(f"  Element prices: {sorted(prices.items())}")


if __name__ == "__main__":
    demo_set_cover()
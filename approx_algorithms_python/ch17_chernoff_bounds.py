"""
Chapter 17: Chernoff Bounds (Set Cover Randomized Rounding)
===========================================================
Vazirani Ch. 17: Chernoff Bounds applied to Set Cover Randomized Rounding.
Implements:
1. Solving Set Cover LP relaxation.
2. Randomized rounding (multiround selection) to ensure high probability cover.
3. Analysis of cover success rate and cost ratio compared to deterministic rounding.
"""

import math
import random
from typing import List, Set, Tuple, Dict
from ch12_14_lp_algorithms import Simplex

def solve_set_cover_lp(
    universe: Set[int],
    sets: Dict[int, Set[int]],  # map set_idx -> set of elements
    costs: Dict[int, float]
) -> Tuple[List[float], float]:
    """
    Solve the Set Cover LP relaxation:
    Minimize sum_s cost_s * x_s
    s.t.
      sum_{s: e in s} x_s >= 1 for all elements e in universe
      x_s >= 0
    """
    elements = sorted(list(universe))
    n_sets = len(sets)
    n_elements = len(elements)
    
    # We formulate standard form: Ax <= b, x >= 0
    # In standard Set Cover, constraints are sum_{s: e in s} x_s >= 1
    # We solve using the Dual LP to ensure feasibility (b_dual >= 0)
    # Primal: Min c^T x s.t. M x >= 1, x >= 0
    # Dual: Max 1^T y s.t. M^T y <= c, y >= 0
    # Here M is the element-set incidence matrix (n_elements x n_sets)
    # So M^T is (n_sets x n_elements)
    
    A_dual = []
    b_dual = []
    c_dual = [1.0] * n_elements
    
    for s_idx in sorted(sets.keys()):
        row = [0.0] * n_elements
        for j, elem in enumerate(elements):
            if elem in sets[s_idx]:
                row[j] = 1.0
        A_dual.append(row)
        b_dual.append(costs[s_idx])
        
    solver = Simplex(A_dual, b_dual, c_dual)
    dual_sol, dual_obj = solver.solve()
    
    if dual_sol is None:
        return [0.0] * n_sets, 0.0
        
    # Extract primal solution from slack variables coefficients in final objective row
    # In final objective row, coefficients of slack variables correspond to primal variables.
    # Col indices of slacks are self.n_cols (n_elements) to n_elements + n_sets
    final_obj = solver.obj_row
    primal_sol = []
    for s_idx in range(n_sets):
        val = final_obj[n_elements + s_idx]
        primal_sol.append(max(0.0, val))
        
    return primal_sol, dual_obj

def set_cover_randomized_rounding(
    universe: Set[int],
    sets: Dict[int, Set[int]],
    costs: Dict[int, float],
    x_lp: List[float],
    c_factor: float = 1.5
) -> Tuple[Set[int], float, bool]:
    """
    Applies randomized rounding to LP Set Cover solution.
    Repeats selection for t = c * ln(n) rounds to ensure high probability cover.
    
    Returns:
        chosen_sets: Indices of selected sets.
        cost: Total cost of selected sets.
        is_valid: True if all elements are covered.
    """
    n_elements = len(universe)
    n_sets = len(sets)
    
    # Number of rounds t = c * ln(n_elements)
    t = int(ceil_val := math.ceil(c_factor * math.log(max(2, n_elements))))
    
    chosen_sets = set()
    
    # Run t independent rounds
    for _ in range(t):
        for s_idx in range(n_sets):
            if random.random() < x_lp[s_idx]:
                chosen_sets.add(s_idx)
                
    # Calculate coverage and cost
    covered = set()
    for s_idx in chosen_sets:
        covered.update(sets[s_idx])
        
    cost = sum(costs[s_idx] for s_idx in chosen_sets)
    is_valid = (covered == universe)
    
    return chosen_sets, cost, is_valid

def demo_chernoff_bounds():
    print("=" * 60)
    print("Chapter 17: Set Cover via Randomized Rounding")
    print("=" * 60)
    
    # Instance: 15 elements, 6 sets
    universe = set(range(15))
    sets = {
        0: {0, 1, 2, 3, 4},
        1: {3, 4, 5, 6, 7},
        2: {6, 7, 8, 9, 10},
        3: {9, 10, 11, 12, 13},
        4: {12, 13, 14, 0, 1},
        5: {2, 5, 8, 11, 14}
    }
    costs = {0: 2.0, 1: 2.0, 2: 2.0, 3: 2.0, 4: 2.0, 5: 1.5}
    
    print("\n1. Set Cover Instance:")
    print(f"  Universe Size: {len(universe)}")
    print("  Available Sets & Costs:")
    for k, v in sets.items():
        print(f"    Set {k}: {sorted(list(v))} (cost: {costs[k]})")
        
    x_lp, lp_obj = solve_set_cover_lp(universe, sets, costs)
    print(f"\n  LP Relaxation Optimal Value: {lp_obj:.4f}")
    print(f"  LP Variable Solution x*:     {[round(xi, 4) for xi in x_lp]}")
    
    # Run randomized rounding trials
    print("\n2. Randomized Rounding Simulation (varying scaling factor c):")
    for c in [0.5, 1.0, 1.5, 2.0]:
        t = math.ceil(c * math.log(len(universe)))
        success_count = 0
        total_cost = 0.0
        n_trials = 200
        
        for _ in range(n_trials):
            _, cost, is_valid = set_cover_randomized_rounding(universe, sets, costs, x_lp, c_factor=c)
            if is_valid:
                success_count += 1
                total_cost += cost
                
        avg_cost = total_cost / max(1, success_count)
        success_rate = success_count / n_trials
        print(f"  c={c:.1f} (rounds={t}): success rate={success_rate:.2%}, avg cost of valid covers={avg_cost:.2f}")

if __name__ == "__main__":
    demo_chernoff_bounds()

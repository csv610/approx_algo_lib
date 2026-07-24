"""
Chapter 16: Randomized Rounding
================================
Vazirani Ch. 16: Randomized Rounding for Max-SAT.
Implements:
1. Max-SAT LP relaxation solver using the Simplex class.
2. Randomized rounding assignment.
3. Coin-flip baseline (0.5 probability).
4. Analysis and comparison of average satisfied weights.
"""

import random
from typing import List, Tuple, Set
from ch12_14_lp_algorithms import Simplex

def solve_max_sat_lp(
    n_vars: int,
    clauses: List[Tuple[Set[int], Set[int]]],  # (positive_lits, negative_lits)
    weights: List[float]
) -> Tuple[List[float], List[float], float]:
    """
    Formulates and solves the LP relaxation for Max-SAT:
    Maximize sum_j w_j * z_j
    s.t.
      z_j <= sum_{i in C_j^+} y_i + sum_{i in C_j^-} (1 - y_i)
      0 <= y_i <= 1
      0 <= z_j <= 1
      
    Returns:
        y: LP values for variables.
        z: LP values for clauses.
        obj: LP optimal objective value.
    """
    n_clauses = len(clauses)
    # Total LP variables: n_vars (y_0..y_{n-1}) + n_clauses (z_0..z_{m-1})
    # Total column count: n_vars + n_clauses
    n_cols = n_vars + n_clauses
    
    A = []
    b = []
    
    # 1. z_j - sum_{i in C_j^+} y_i + sum_{i in C_j^-} y_i <= |C_j^-|
    for j, (pos, neg) in enumerate(clauses):
        row = [0.0] * n_cols
        # z_j coefficient
        row[n_vars + j] = 1.0
        # y_i coefficients
        for i in pos:
            row[i] = -1.0
        for i in neg:
            row[i] = 1.0
            
        A.append(row)
        b.append(float(len(neg)))
        
    # 2. y_i <= 1
    for i in range(n_vars):
        row = [0.0] * n_cols
        row[i] = 1.0
        A.append(row)
        b.append(1.0)
        
    # 3. z_j <= 1
    for j in range(n_clauses):
        row = [0.0] * n_cols
        row[n_vars + j] = 1.0
        A.append(row)
        b.append(1.0)
        
    # Objective: Maximize sum_j w_j * z_j
    c = [0.0] * n_cols
    for j in range(n_clauses):
        c[n_vars + j] = weights[j]
        
    # Solve
    solver = Simplex(A, b, c)
    sol, obj = solver.solve()
    
    if sol is None:
        return [0.0] * n_vars, [0.0] * n_clauses, 0.0
        
    y = sol[:n_vars]
    z = sol[n_vars:]
    return y, z, obj

def evaluate_assignment(
    assignment: List[bool],
    clauses: List[Tuple[Set[int], Set[int]]],
    weights: List[float]
) -> float:
    """Calculate the total weight of satisfied clauses."""
    total_weight = 0.0
    for j, (pos, neg) in enumerate(clauses):
        satisfied = False
        for i in pos:
            if assignment[i]:
                satisfied = True
                break
        if not satisfied:
            for i in neg:
                if not assignment[i]:
                    satisfied = True
                    break
        if satisfied:
            total_weight += weights[j]
    return total_weight

def randomized_rounding_max_sat(
    n_vars: int,
    clauses: List[Tuple[Set[int], Set[int]]],
    weights: List[float],
    y_lp: List[float],
    trials: int = 200
) -> Tuple[List[bool], float, float]:
    """
    Applies randomized rounding to LP variable values.
    
    Returns:
        best_assignment: The assignment that achieved the highest weight.
        avg_weight: The average weight over all trials.
        best_weight: The maximum weight achieved.
    """
    best_weight = -1.0
    best_assignment = []
    total_weight = 0.0
    
    for _ in range(trials):
        assignment = []
        for i in range(n_vars):
            assignment.append(random.random() < y_lp[i])
            
        weight = evaluate_assignment(assignment, clauses, weights)
        total_weight += weight
        if weight > best_weight:
            best_weight = weight
            best_assignment = assignment
            
    return best_assignment, total_weight / trials, best_weight

def coin_flip_max_sat(
    n_vars: int,
    clauses: List[Tuple[Set[int], Set[int]]],
    weights: List[float],
    trials: int = 200
) -> Tuple[List[bool], float, float]:
    """Random coin-flip baseline (set variable to True with prob 0.5)."""
    best_weight = -1.0
    best_assignment = []
    total_weight = 0.0
    
    for _ in range(trials):
        assignment = [random.random() < 0.5 for _ in range(n_vars)]
        weight = evaluate_assignment(assignment, clauses, weights)
        total_weight += weight
        if weight > best_weight:
            best_weight = weight
            best_assignment = assignment
            
    return best_assignment, total_weight / trials, best_weight

def demo_randomized_rounding():
    print("=" * 60)
    print("Chapter 16: Randomized Rounding for Max-SAT")
    print("=" * 60)
    
    # 3 variables, 4 clauses
    # C_0 = x_0 or x_1
    # C_1 = not x_1 or x_2
    # C_2 = not x_0 or not x_2
    # C_3 = x_1
    n_vars = 3
    clauses = [
        ({0, 1}, set()),      # x_0 or x_1
        ({2}, {1}),          # not x_1 or x_2
        (set(), {0, 2}),      # not x_0 or not x_2
        ({1}, set())         # x_1
    ]
    weights = [1.0, 2.0, 1.5, 3.0]
    
    print("\n1. Max-SAT CNF Instance:")
    print("  Variables: x_0, x_1, x_2")
    print("  Clauses & Weights:")
    for j, (pos, neg) in enumerate(clauses):
        lits = [f"x_{i}" for i in pos] + [f"not x_{i}" for i in neg]
        print(f"    C_{j}: ({' or '.join(lits)}), weight: {weights[j]}")
        
    y_lp, z_lp, lp_obj = solve_max_sat_lp(n_vars, clauses, weights)
    print(f"\n  LP Optimal Objective Value: {lp_obj:.4f}")
    print(f"  LP Variable Solution y*:    {[round(yi, 4) for yi in y_lp]}")
    print(f"  LP Clause Solution z*:      {[round(zj, 4) for zj in z_lp]}")
    
    # Run algorithms
    best_rr, avg_rr, max_rr = randomized_rounding_max_sat(n_vars, clauses, weights, y_lp, trials=500)
    best_cf, avg_cf, max_cf = coin_flip_max_sat(n_vars, clauses, weights, trials=500)
    
    # Find theoretical best by checking all 8 assignments
    opt_val = 0.0
    best_all = None
    for bitmask in range(8):
        assignment = [(bitmask & (1 << i)) > 0 for i in range(3)]
        val = evaluate_assignment(assignment, clauses, weights)
        if val > opt_val:
            opt_val = val
            best_all = assignment
            
    print(f"\n  Exact Optimal Value:        {opt_val:.2f} (assignment: {best_all})")
    print(f"  Randomized Rounding (avg):  {avg_rr:.4f} (ratio to OPT: {avg_rr/opt_val:.4f}, bound: 1-1/e ~ 0.632)")
    print(f"  Randomized Rounding (best): {max_rr:.4f}")
    print(f"  Coin Flip Baseline (avg):   {avg_cf:.4f} (ratio to OPT: {avg_cf/opt_val:.4f}, bound: 0.5)")
    print(f"  Coin Flip Baseline (best):  {max_cf:.4f}")

if __name__ == "__main__":
    demo_randomized_rounding()

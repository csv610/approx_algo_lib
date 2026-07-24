"""
Chapter 26: Semidefinite Programming (Max 2-SAT)
=================================================
Vazirani Ch. 26: SDP-based 0.878-approximation for Max 2-SAT.
Implements:
1. Vector representation for variables and a special "True" vector v_0.
2. Spherical gradient projection optimizer for Max 2-SAT SDP relaxation.
3. Hyperplane rounding with sign matching against v_0.
4. Performance analysis on test instances.
"""

import math
import random
from typing import List, Tuple, Dict, Set

def evaluate_sat_assignment(
    assignment: List[bool],  # 1-indexed, size n_vars + 1
    clauses: List[Tuple[int, int]]  # literals: positive if >0, negative if <0
) -> List[bool]:
    """Returns a list of booleans indicating if each clause is satisfied."""
    satisfied = []
    for l1, l2 in clauses:
        sat1 = assignment[abs(l1)] if l1 > 0 else not assignment[abs(l1)]
        if l2 == 0:  # single literal clause
            satisfied.append(sat1)
        else:
            sat2 = assignment[abs(l2)] if l2 > 0 else not assignment[abs(l2)]
            satisfied.append(sat1 or sat2)
    return satisfied

def optimize_max_2sat_vectors(
    n_vars: int,
    clauses: List[Tuple[int, int]],
    weights: List[float],
    dim: int = 8,
    lr: float = 0.05,
    epochs: int = 300
) -> List[List[float]]:
    """
    Optimizes vectors v_0, v_1, ..., v_{n_vars} on the unit sphere in R^dim.
    v_0 corresponds to the "True" anchor vector.
    
    We maximize the SDP relaxation:
    Maximize sum_{C_j} w_j * val(C_j)
    where:
      - For C_j = (l1 or l2):
        val = 1/4 * (3 + sgn(l1)*v_{|l1|} . v_0 + sgn(l2)*v_{|l2|} . v_0 - sgn(l1)*sgn(l2)*v_{|l1|} . v_{|l2|})
      - For C_j = (l1):
        val = 1/2 * (1 + sgn(l1)*v_{|l1|} . v_0)
    """
    n_vectors = n_vars + 1 # v_0 to v_n
    
    # Initialize random unit vectors
    v = []
    for _ in range(n_vectors):
        vec = [random.gauss(0.0, 1.0) for _ in range(dim)]
        mag = math.sqrt(sum(x*x for x in vec))
        v.append([x/mag for x in vec])
        
    for _ in range(epochs):
        new_v = []
        for i in range(n_vectors):
            grad = [0.0] * dim
            # Compute gradient dF/dv_i
            for j, (l1, l2) in enumerate(clauses):
                w = weights[j]
                s1 = 1 if l1 > 0 else -1
                v1_idx = abs(l1)
                
                if l2 == 0:  # Single literal: 1/2 * (1 + s1 * v_{v1_idx} . v_0)
                    if i == 0: # v_0
                        for k in range(dim):
                            grad[k] += 0.5 * w * s1 * v[v1_idx][k]
                    elif i == v1_idx: # v_i
                        for k in range(dim):
                            grad[k] += 0.5 * w * s1 * v[0][k]
                else:  # Two literals
                    s2 = 1 if l2 > 0 else -1
                    v2_idx = abs(l2)
                    
                    if i == 0: # v_0
                        for k in range(dim):
                            grad[k] += 0.25 * w * (s1 * v[v1_idx][k] + s2 * v[v2_idx][k])
                    elif i == v1_idx: # v_{|l1|}
                        for k in range(dim):
                            grad[k] += 0.25 * w * (s1 * v[0][k] - s1 * s2 * v[v2_idx][k])
                    elif i == v2_idx: # v_{|l2|}
                        for k in range(dim):
                            grad[k] += 0.25 * w * (s2 * v[0][k] - s1 * s2 * v[v1_idx][k])
                            
            # Update: Gradient Ascent (since we are maximizing!)
            updated = [v[i][k] + lr * grad[k] for k in range(dim)]
            mag = math.sqrt(sum(x*x for x in updated))
            if mag < 1e-9:
                new_v.append(v[i])
            else:
                new_v.append([x/mag for x in updated])
        v = new_v
    return v

def goemans_williamson_max_2sat(
    n_vars: int,
    clauses: List[Tuple[int, int]],
    weights: List[float],
    vectors: List[List[float]],
    trials: int = 500
) -> Tuple[List[bool], float, float]:
    """
    Applies hyperplane rounding for Max 2-SAT.
    Generates a random vector r.
    Variable x_i is set to True if sign(v_i . r) == sign(v_0 . r), else False.
    
    Returns:
        best_assignment: boolean assignment list (1-indexed, size n_vars+1).
        avg_weight: average weight satisfied.
        best_weight: maximum weight satisfied.
    """
    dim = len(vectors[0])
    best_weight = -1.0
    best_assignment = []
    total_weight = 0.0
    
    for _ in range(trials):
        r = [random.gauss(0.0, 1.0) for _ in range(dim)]
        mag = math.sqrt(sum(x*x for x in r))
        r = [x/mag for x in r]
        
        # Determine signs of dot products
        dot_v0 = sum(vectors[0][k] * r[k] for k in range(dim))
        sign_v0 = dot_v0 >= 0.0
        
        assignment = [False] * (n_vars + 1)
        for i in range(1, n_vars + 1):
            dot_vi = sum(vectors[i][k] * r[k] for k in range(dim))
            sign_vi = dot_vi >= 0.0
            assignment[i] = (sign_vi == sign_v0)
            
        sat_mask = evaluate_sat_assignment(assignment, clauses)
        weight = sum(w for sat, w in zip(sat_mask, weights) if sat)
        total_weight += weight
        
        if weight > best_weight:
            best_weight = weight
            best_assignment = assignment
            
    return best_assignment, total_weight / trials, best_weight

def demo_max_2sat():
    print("=" * 60)
    print("Chapter 26: SDP for Maximum 2-SAT")
    print("=" * 60)
    
    # 2 variables, 4 clauses (contradictory set)
    # C_0 = x_1 or x_2
    # C_1 = not x_1 or x_2
    # C_2 = x_1 or not x_2
    # C_3 = not x_1 or not x_2
    n_vars = 2
    clauses = [(1, 2), (-1, 2), (1, -2), (-1, -2)]
    weights = [1.0, 1.0, 1.0, 1.0]
    
    print("\n1. 2-Variable Max 2-SAT (All Combinations):")
    print("  Clauses & Weights:")
    for j, (l1, l2) in enumerate(clauses):
        lit1 = f"x_{l1}" if l1 > 0 else f"not x_{-l1}"
        lit2 = f"x_{l2}" if l2 > 0 else f"not x_{-l2}"
        print(f"    C_{j}: ({lit1} or {lit2}), weight: {weights[j]}")
        
    vectors = optimize_max_2sat_vectors(n_vars, clauses, weights, dim=3, lr=0.1, epochs=300)
    print("\n  Optimized Vector Embeddings:")
    print(f"    v_0 (True Anchor): {[round(x, 4) for x in vectors[0]]}")
    for i in range(1, n_vars + 1):
        print(f"    v_{i} (Var x_{i}):   {[round(x, 4) for x in vectors[i]]}")
        
    best_assign, avg_w, max_w = goemans_williamson_max_2sat(n_vars, clauses, weights, vectors, trials=1000)
    print(f"\n  Exact Optimal Max 2-SAT Weight: 3.0")
    print(f"  GW Rounding Average Weight:    {avg_w:.4f} (ratio to OPT: {avg_w/3.0:.4f}, bound: 0.878)")
    print(f"  GW Rounding Best Weight:       {max_w:.4f} (assignment: {best_assign[1:]})")
    
    # Example 2: Larger Instance (5 variables, 8 clauses)
    n_vars2 = 5
    clauses2 = [
        (1, 2), (-2, 3), (-3, 4), (4, 5), (-5, -1),
        (1, 0), (-2, 0), (3, -5)
    ]
    weights2 = [1.5, 2.0, 1.0, 1.5, 2.0, 3.0, 1.0, 2.5]
    
    print("\n2. Larger Max 2-SAT Instance (5 variables, 8 clauses):")
    vectors2 = optimize_max_2sat_vectors(n_vars2, clauses2, weights2, dim=8, lr=0.05, epochs=300)
    best_assign2, avg_w2, max_w2 = goemans_williamson_max_2sat(n_vars2, clauses2, weights2, vectors2, trials=1000)
    
    # Solve exactly by brute forcing 32 combinations
    opt_w2 = 0.0
    for bitmask in range(32):
        assignment = [False] + [(bitmask & (1 << i)) > 0 for i in range(5)]
        sat_mask = evaluate_sat_assignment(assignment, clauses2)
        w = sum(wt for sat, wt in zip(sat_mask, weights2) if sat)
        if w > opt_w2:
            opt_w2 = w
            
    print(f"  Exact Optimal Max 2-SAT Weight: {opt_w2:.2f}")
    print(f"  GW Rounding Average Weight:    {avg_w2:.4f} (ratio to OPT: {avg_w2/opt_w2:.4f}, bound: 0.878)")
    print(f"  GW Rounding Best Weight:       {max_w2:.4f} (assignment: {best_assign2[1:]})")

if __name__ == "__main__":
    demo_max_2sat()

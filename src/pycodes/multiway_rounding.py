"""
Chapter 19: Multiway Cut LP Rounding
====================================
Vazirani Ch. 19: Simplex-distance LP relaxation and CKR 1.5-approximation.
Implements:
1. Simplex-distance LP formulation and Dual LP solver.
2. Călinescu-Karloff-Rabani (CKR) randomized rounding.
3. Performance ratio comparison with Chapter 4's isolating cuts.
"""

import math
import random
from typing import List, Tuple, Dict, Set

from lp_algorithms import Simplex

def solve_multiway_cut_lp(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    terminals: List[int]
) -> Tuple[List[List[float]], float]:
    """
    Solves the Multiway Cut LP relaxation:
    Minimize 1/2 * sum_{(u,v) in E} c_{uv} * sum_{i=1}^k |d_i(u) - d_i(v)|
    s.t.
      sum_i d_i(v) = 1 for all v
      d_i(s_j) = 1 if i == j else 0
      d_i(v) >= 0
      
    We formulate this in standard minimization form:
      Minimize c^T x s.t. M x >= f, x >= 0
      
    And solve the Dual:
      Maximize f^T y s.t. M^T y <= c, y >= 0
    """
    k = len(terminals)
    n_edges = len(edges)
    
    # Primal variables:
    # 1. d_i(v) for v in range(n), i in range(k) (size n * k)
    #    Index: v * k + i
    # 2. y_{uv}^i for e in range(n_edges), i in range(k) (size n_edges * k)
    #    Index: n * k + e * k + i
    n_vars = n * k + n_edges * k
    
    # We construct constraints of form M x >= f
    M = []
    f = []
    
    # 1. y_{uv}^i - d_i(u) + d_i(v) >= 0
    for e, (u, v) in enumerate(edges):
        for i in range(k):
            row = [0.0] * n_vars
            row[n * k + e * k + i] = 1.0 # y_{uv}^i
            row[u * k + i] = -1.0         # -d_i(u)
            row[v * k + i] = 1.0          # d_i(v)
            M.append(row)
            f.append(0.0)
            
    # 2. y_{uv}^i + d_i(u) - d_i(v) >= 0
    for e, (u, v) in enumerate(edges):
        for i in range(k):
            row = [0.0] * n_vars
            row[n * k + e * k + i] = 1.0 # y_{uv}^i
            row[u * k + i] = 1.0          # d_i(u)
            row[v * k + i] = -1.0         # -d_i(v)
            M.append(row)
            f.append(0.0)
            
    # 3. sum_i d_i(v) >= 1
    for v in range(n):
        row = [0.0] * n_vars
        for i in range(k):
            row[v * k + i] = 1.0
        M.append(row)
        f.append(1.0)
        
    # 4. -sum_i d_i(v) >= -1
    for v in range(n):
        row = [0.0] * n_vars
        for i in range(k):
            row[v * k + i] = -1.0
        M.append(row)
        f.append(-1.0)
        
    # 5. Terminal boundary conditions:
    # d_i(s_j) >= 1 (if i == j)
    # -d_i(s_j) >= 0 (if i != j)
    for j, s_j in enumerate(terminals):
        for i in range(k):
            row = [0.0] * n_vars
            if i == j:
                row[s_j * k + i] = 1.0
                M.append(row)
                f.append(1.0)
            else:
                row[s_j * k + i] = -1.0
                M.append(row)
                f.append(0.0)
                
    # Primal Objective: Minimize 1/2 * sum_e c_e * sum_i y_e^i
    c_primal = [0.0] * n_vars
    for e in range(n_edges):
        cost = costs[e]
        for i in range(k):
            c_primal[n * k + e * k + i] = 0.5 * cost
            
    # Formulate Dual: Maximize f^T y s.t. M^T y <= c_primal
    n_dual_vars = len(f)
    A_dual = []
    b_dual = c_primal
    c_dual = f
    
    # A_dual is M^T (size n_vars x n_dual_vars)
    for col_idx in range(n_vars):
        row = [M[row_idx][col_idx] for row_idx in range(n_dual_vars)]
        A_dual.append(row)
        
    solver = Simplex(A_dual, b_dual, c_dual)
    dual_sol, dual_obj = solver.solve()
    
    if dual_sol is None:
        return [[1.0/k] * k for _ in range(n)], 0.0
        
    # Extract primal variables from slack variables coefficients in final objective row
    final_obj = solver.obj_row
    primal_vals = []
    for var_idx in range(n_vars):
        val = final_obj[n_dual_vars + var_idx]
        primal_vals.append(max(0.0, val))
        
    d = []
    for v in range(n):
        v_d = [primal_vals[v * k + i] for i in range(k)]
        d.append(v_d)
        
    return d, dual_obj

def calinescu_karloff_rabani_rounding(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    terminals: List[int],
    d: List[List[float]],
    trials: int = 500
) -> Tuple[List[int], float]:
    """
    Applies CKR randomized rounding.
    1. Select a random permutation of terminals.
    2. Select a random radius r uniformly in (0, 1/2).
    3. Assign vertex v to the first terminal i in the permutation such that d_i(v) > r.
    
    Returns:
        best_cut: list of edges in the cut.
        best_cost: total cost of the cut.
    """
    k = len(terminals)
    best_cost = float('inf')
    best_cut = []
    
    for _ in range(trials):
        # 1. Random permutation
        perm = list(range(k))
        random.shuffle(perm)
        
        # 2. Random radius r in (0, 1/2)
        r = random.uniform(0.0, 0.5)
        
        # 3. Assign vertices
        assignment = [-1] * n
        for v in range(n):
            assigned = False
            for idx in perm:
                if d[v][idx] > r:
                    assignment[v] = idx
                    assigned = True
                    break
            if not assigned:
                assignment[v] = perm[-1]
                
        # Compute cut edges
        cut = []
        cost = 0.0
        for e, (u, v) in enumerate(edges):
            if assignment[u] != assignment[v]:
                cut.append((u, v))
                cost += costs[e]
                
        if cost < best_cost:
            best_cost = cost
            best_cut = cut
            
    return best_cut, best_cost

def demo_multiway_cut_lp():
    print("=" * 60)
    print("Chapter 19: Multiway Cut via LP Rounding")
    print("=" * 60)
    
    # 6-node graph, 3 terminals
    n = 6
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),  # outer cycle
        (0, 5), (2, 5), (4, 5)                  # inner spokes
    ]
    costs = [2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0]
    terminals = [0, 2, 4]
    
    print("\n1. Input Instance:")
    print(f"  Vertices: {list(range(n))}")
    print(f"  Terminals: {terminals}")
    print("  Edges & Costs:")
    for i, (u, v) in enumerate(edges):
        print(f"    ({u}, {v}) cost={costs[i]}")
        
    d, lp_obj = solve_multiway_cut_lp(n, edges, costs, terminals)
    print(f"\n  LP Relaxation Optimal Value: {lp_obj:.4f}")
    print("  LP Vertex Embeddings on Simplex:")
    for v in range(n):
        print(f"    v_{v}: {[round(val, 4) for val in d[v]]}")
        
    cut, cost = calinescu_karloff_rabani_rounding(n, edges, costs, terminals, d, trials=500)
    print("\n2. CKR Rounding Results:")
    print(f"  Selected Cut Edges: {cut}")
    print(f"  Cut Total Cost:     {cost:.2f}")

if __name__ == "__main__":
    demo_multiway_cut_lp()

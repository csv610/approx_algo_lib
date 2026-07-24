"""
Chapter 22: Steiner Network (Kamal Jain's Iterative Rounding)
============================================================
Vazirani Ch. 23: Jain's Iterative Rounding 2-approximation for SNDP.
Implements:
1. Formulation of the Dual LP for Survivable Network Design.
2. Iterative LP solving and variable rounding.
3. Updated connectivity requirements at each phase.
"""

import math
from typing import List, Tuple, Dict, Set

from ch12_14_lp_algorithms import Simplex

def get_all_cuts(n: int) -> List[Set[int]]:
    """Generates all non-empty proper cuts in a graph of size n."""
    cuts = []
    # Loop over all subsets of vertices (excluding empty and full)
    for i in range(1, 1 << (n - 1)):
        cut = set()
        for j in range(n):
            if (i & (1 << j)) > 0:
                cut.add(j)
        cuts.append(cut)
    return cuts

def solve_sndp_lp_phase(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    r: Dict[Tuple[int, int], int],
    fixed_edges: Set[int]  # Indices of edges rounded to 1
) -> List[float]:
    """
    Solves the continuous LP relaxation for SNDP under current requirements.
    We formulate the Dual LP:
      Maximize sum_S f(S) * y_S - sum_e z_e
      s.t.
        sum_{S: e in delta(S)} y_S - z_e <= c_e  (for all e not in fixed_edges)
        y_S >= 0, z_e >= 0
        
    Note: For e in fixed_edges, x_e is fixed to 1, so it does not participate.
    The connectivity requirement f(S) is decreased by the number of fixed edges crossing delta(S).
    """
    cuts = get_all_cuts(n)
    n_cuts = len(cuts)
    n_edges = len(edges)
    
    # Compute active edge indices
    active_edges = [i for i in range(n_edges) if i not in fixed_edges]
    n_active = len(active_edges)
    
    # Calculate requirements f(S) for each cut
    f = []
    for cut in cuts:
        # Initial requirement: max_{u in S, v not in S} r(u,v)
        val = 0
        for (u, v), req in r.items():
            u_in = u in cut
            v_in = v in cut
            if u_in != v_in:
                val = max(val, req)
                
        # Subtract fixed edges crossing this cut
        fixed_crossing = 0
        for edge_idx in fixed_edges:
            u, v = edges[edge_idx]
            if (u in cut) != (v in cut):
                fixed_crossing += 1
                
        f.append(max(0, val - fixed_crossing))
        
    # Dual LP Variables: y_S for S in cuts (index 0..n_cuts-1) and z_e for e in active_edges (index n_cuts..n_cuts+n_active-1)
    n_cols = n_cuts + n_active
    
    A_dual = []
    b_dual = []
    
    # Constraints: for each active edge e, sum_{S: e in delta(S)} y_S - z_e <= c_e
    for j, edge_idx in enumerate(active_edges):
        u, v = edges[edge_idx]
        row = [0.0] * n_cols
        # y_S coefficients
        for s_idx, cut in enumerate(cuts):
            if (u in cut) != (v in cut):
                row[s_idx] = 1.0
        # z_e coefficient
        row[n_cuts + j] = -1.0
        
        A_dual.append(row)
        b_dual.append(costs[edge_idx])
        
    # Objective: Maximize sum_S f(S) * y_S - sum_e 1 * z_e
    c_dual = [float(val) for val in f] + [-1.0] * n_active
    
    solver = Simplex(A_dual, b_dual, c_dual)
    dual_sol, dual_obj = solver.solve()
    
    if dual_sol is None:
        return [0.0] * n_edges
        
    # Extract primal variables x_e from slack variables coefficients in final objective row
    final_obj = solver.obj_row
    x = [0.0] * n_edges
    
    # Fixed edges are 1.0
    for idx in fixed_edges:
        x[idx] = 1.0
        
    # Active edges primal values are final_obj[n_cols + j]
    for j, edge_idx in enumerate(active_edges):
        val = final_obj[n_cols + j]
        x[edge_idx] = max(0.0, val)
        
    return x

def jain_iterative_rounding(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    r: Dict[Tuple[int, int], int]
) -> List[Tuple[int, int]]:
    """
    Kamal Jain's Iterative Rounding Algorithm.
    Repeatedly solves LP, rounds any edge with x_e >= 0.5 to 1, and re-solves.
    """
    fixed_edges = set()
    
    while True:
        # Solve LP phase
        x = solve_sndp_lp_phase(n, edges, costs, r, fixed_edges)
        
        # Check if we are done (all requirements met)
        cuts = get_all_cuts(n)
        requirements_satisfied = True
        for cut in cuts:
            # Check requirement f(S)
            req = 0
            for (u, v), val in r.items():
                if (u in cut) != (v in cut):
                    req = max(req, val)
            if req == 0:
                continue
            # Check actual capacity in current fixed sets
            capacity = sum(1 for idx in fixed_edges if (edges[idx][0] in cut) != (edges[idx][1] in cut))
            if capacity < req:
                requirements_satisfied = False
                break
                
        if requirements_satisfied:
            break
            
        # Find active edge with highest x_e
        best_val = -1.0
        best_idx = -1
        for idx in range(len(edges)):
            if idx not in fixed_edges:
                if x[idx] > best_val:
                    best_val = x[idx]
                    best_idx = idx
                    
        # Jain's theorem guarantees there exists an edge with x_e >= 0.5
        # Round it to 1
        if best_idx != -1 and best_val >= 0.4999:
            fixed_edges.add(best_idx)
        else:
            # Fallback if numerical issues occur
            if best_idx != -1:
                fixed_edges.add(best_idx)
            else:
                break
                
    return [edges[idx] for idx in fixed_edges]

def demo_steiner_network():
    print("=" * 60)
    print("Chapter 22: Steiner Network (Jain's Iterative Rounding)")
    print("=" * 60)
    
    # 4-node graph
    n = 4
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
    costs = [2.0, 3.0, 2.0, 3.0, 4.0]
    
    # Connectivity requirements:
    # r(0,2) = 2 (need 2 disjoint paths)
    # r(1,3) = 1 (need 1 path)
    r = {
        (0, 2): 2,
        (1, 3): 1
    }
    
    print("\n1. Input Graph and Requirements:")
    print("  Edges & Costs:")
    for i, (u, v) in enumerate(edges):
        print(f"    Edge {i}: ({u}, {v}) cost={costs[i]}")
    print("  Connectivity Requirements:")
    for (u, v), req in r.items():
        print(f"    r({u}, {v}) = {req}")
        
    chosen = jain_iterative_rounding(n, edges, costs, r)
    print("\n2. Iterative Rounding Result:")
    print(f"  Selected Edges: {chosen}")
    total_cost = sum(costs[edges.index(e)] if e in edges else costs[edges.index((e[1], e[0]))] for e in chosen)
    print(f"  Total Network Cost: {total_cost:.2f}")

if __name__ == "__main__":
    demo_steiner_network()

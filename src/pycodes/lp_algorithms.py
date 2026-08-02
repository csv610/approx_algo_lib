"""
Chapter 12-14: LP-Duality Based Algorithms
===========================================
Vazirani Part II: LP-based approximation algorithms.

Key techniques:
- Rounding (Ch. 14): Solve LP relaxation, round fractional solution
- Primal-Dual Schema (Ch. 13, 22): Construct primal/dual solutions simultaneously
- Dual Fitting (Ch. 13): Analyze greedy via dual fitting

This chapter implements:
- Set Cover via LP Rounding (Ch. 14)
- Set Cover via Primal-Dual (Ch. 13/22)
- Vertex Cover via LP Rounding (2-approx)
"""

from typing import Dict, List, Set, Tuple, Optional
import itertools
import math
class Simplex:
    """Basic simplex implementation for small LPs.
    
    Maximize c^T x subject to Ax <= b, x >= 0
    A: m x n matrix
    b: m vector
    c: n vector
    """
    
    def __init__(self, A: List[List[float]], b: List[float], c: List[float]):
        self.m = len(A)
        self.n = len(A[0]) if A else 0
        self.A = [row[:] for row in A]
        self.b = b[:]
        self.c = c[:]
        self.obj_row = []
    
    def solve(self) -> Tuple[List[float], float]:
        """Returns (solution, optimal_value) or (None, None) if unbounded/infeasible."""
        # Add slack variables for standard form
        # Maximize c^T x s.t. Ax + s = b, x,s >= 0
        
        # Initial tableau
        # Rows: constraints + objective
        # Cols: x1..xn, s1..sm, RHS
        tab = []
        for i in range(self.m):
            row = self.A[i] + [1 if j == i else 0 for j in range(self.m)] + [self.b[i]]
            tab.append(row)
        obj = [-ci for ci in self.c] + [0] * self.m + [0]
        tab.append(obj)
        
        # Simplex iterations
        while True:
            # Find entering variable (most negative in obj row)
            entering = -1
            min_val = -1e-9
            for j in range(self.n + self.m):
                if tab[-1][j] < min_val:
                    min_val = tab[-1][j]
                    entering = j
            if entering == -1:
                break  # Optimal
            
            # Find leaving variable (min ratio test)
            leaving = -1
            min_ratio = float('inf')
            for i in range(self.m):
                if tab[i][entering] > 1e-9:
                    ratio = tab[i][-1] / tab[i][entering]
                    if ratio < min_ratio:
                        min_ratio = ratio
                        leaving = i
            if leaving == -1:
                return None, float('inf')  # Unbounded
            
            # Pivot
            pivot = tab[leaving][entering]
            tab[leaving] = [v / pivot for v in tab[leaving]]
            
            for i in range(self.m + 1):
                if i != leaving and abs(tab[i][entering]) > 1e-9:
                    factor = tab[i][entering]
                    tab[i] = [tab[i][j] - factor * tab[leaving][j] 
                              for j in range(len(tab[0]))]
        
        # Extract solution
        x = [0.0] * self.n
        for i in range(self.m):
            # Find basic variable in row i: look for column that is
            # approximately a unit vector (1 in row i, ~0 elsewhere)
            best_j = -1
            best_score = -1.0
            for j in range(self.n + self.m):
                if abs(tab[i][j] - 1.0) < 1e-6:
                    # Count how "unit-vector-like" this column is
                    col_sum = sum(abs(tab[r][j]) for r in range(self.m + 1))
                    score = col_sum  # should be ~1.0 for a unit column
                    if score > best_score and score < 1.0 + 1e-4:
                        best_score = score
                        best_j = j
            if best_j >= 0 and best_j < self.n:
                x[best_j] = max(0.0, tab[i][-1])
        
        self.obj_row = tab[-1]
        return x, tab[-1][-1]


# ============================================================
# SET COVER LP RELAXATION
# ============================================================

def set_cover_lp(universe: Set[int], sets: Dict[int, Set[int]], 
                 costs: Dict[int, float]) -> Tuple[List[float], float]:
    """
    Solve LP Relaxation for Set Cover via its DUAL formulation:
    Dual LP: Maximize sum y_e
    Subject to: sum_{e in S} y_e <= c_s for all S
                y_e >= 0
    
    This formulation has b_dual = costs >= 0, which is feasible from the origin y=0.
    The primal variables x_s* are extracted from the coefficients of the slack variables 
    in the final objective row of the dual simplex tableau.
    """
    U = list(universe)
    S = list(sets.keys())
    m = len(S)
    n = len(U)
    
    if m == 0:
        return [0.0] * m, float('inf')
    
    # Dual LP constraints: A_dual y <= b_dual
    A_dual = []
    for s in S:
        row = [1.0 if e in sets[s] else 0.0 for e in U]
        A_dual.append(row)
        
    b_dual = [costs[s] for s in S]
    c_dual = [1.0] * n
    
    simplex = Simplex(A_dual, b_dual, c_dual)
    y, opt = simplex.solve()
    
    if y is None:
        return [0.0] * m, float('inf')
        
    # Primal variables x correspond to slack variables of the dual.
    x = simplex.obj_row[n : n + m]
    return x, opt


def set_cover_lp_rounding(universe: Set[int], sets: Dict[int, Set[int]], 
                          costs: Dict[int, float]) -> Tuple[List[int], float]:
    """
    Set Cover via LP Rounding (Theorem 14.3 in Vazirani).
    
    Algorithm:
    1. Solve LP relaxation
    2. Pick all sets with x_s >= 1/f where f = max frequency
    
    Approximation factor: f (frequency of most frequent element)
    """
    x, opt = set_cover_lp(universe, sets, costs)
    if not x:
        return [], 0.0
    
    S = list(sets.keys())
    m = len(S)
    
    # Compute frequency
    freq = 0
    for e in universe:
        cnt = sum(1 for s in S if e in sets[s])
        freq = max(freq, cnt)
    
    if freq == 0:
        return [], 0.0
    
    threshold = 1.0 / freq
    picked = []
    total_cost = 0.0
    
    for j, s in enumerate(S):
        if x[j] >= threshold:
            picked.append(s)
            total_cost += costs[s]
    
    return picked, total_cost


def set_cover_primal_dual(universe: Set[int], sets: Dict[int, Set[int]], 
                          costs: Dict[int, float]) -> Tuple[List[int], float]:
    """
    Set Cover via Primal-Dual Schema (Algorithm 13.1 / 22.3 in Vazirani).
    
    Primal: min sum c_s x_s
    Dual: max sum y_e
          s.t. sum_{e in S} y_e <= c_s for all S
               y_e >= 0
    
    Algorithm:
    1. y = 0
    2. While there is uncovered element e:
       - Increase y_e until some set S becomes tight (sum_{e in S} y_e = c_s)
       - Add S to cover
       - Mark all elements in S as covered
    
    Approximation factor: f (frequency)
    """
    U = set(universe)
    covered = set()
    y = {e: 0.0 for e in U}
    picked = []
    total_cost = 0.0
    
    # Precompute sets containing each element
    sets_containing = {e: [] for e in U}
    for s, se in sets.items():
        for e in se:
            if e in U:
                sets_containing[e].append(s)
    
    while covered != U:
        # Find uncovered element
        e = next(iter(U - covered))
        
        # Increase y_e until some set becomes tight
        while True:
            # Find minimum slack among sets containing e
            min_slack = float('inf')
            tight_set = None
            
            for s in sets_containing[e]:
                slack = costs[s] - sum(y[ee] for ee in sets[s] if ee in U)
                if slack < min_slack:
                    min_slack = slack
                    tight_set = s
            
            if min_slack <= 1e-9:
                # Already tight (or will be with tiny increase)
                break
            
            # Increase y_e by min_slack
            y[e] += min_slack
            
            # Check if any set became tight
            for s in sets_containing[e]:
                slack = costs[s] - sum(y[ee] for ee in sets[s] if ee in U)
                if abs(slack) < 1e-9:
                    tight_set = s
                    break
            
            if tight_set is not None:
                break
        
        if tight_set is None:
            # Should not happen
            break
        
        # Add tight set to cover
        if tight_set not in picked:
            picked.append(tight_set)
            total_cost += costs[tight_set]
            covered.update(sets[tight_set])
    
    return picked, total_cost


def vertex_cover_lp_rounding(graph: Dict[int, Dict[int, float]]) -> Tuple[Set[int], float]:
    """
    Vertex Cover via LP Rounding (2-approx).
    
    LP: min sum w_v x_v
        s.t. x_u + x_v >= 1 for all (u,v) in E
             x_v >= 0
    
    Rounding: pick all vertices with x_v >= 1/2
    """
    # Map to Set Cover
    edges = []
    for u in graph:
        for v in graph[u]:
            if u < v:
                edges.append((u, v))
                
    universe = set(range(len(edges)))
    edge_to_idx = {e: i for i, e in enumerate(edges)}
    
    vertices = list(graph.keys())
    sets = {}
    for u in vertices:
        sets[u] = set()
        for v in graph[u]:
            e = (min(u, v), max(u, v))
            sets[u].add(edge_to_idx[e])
            
    costs = {v: 1.0 for v in vertices}
    
    x, opt = set_cover_lp(universe, sets, costs)
    
    if opt == float('inf'):
        return set(vertices), float('inf')
        
    cover = {vertices[i] for i in range(len(vertices)) if x[i] >= 0.5}
    return cover, float(len(cover))


def demo_lp_algorithms():
    print("=" * 60)
    print("Chapters 12-14: LP-Duality Based Algorithms")
    print("=" * 60)
    
    # Set Cover example
    print("\n1. Set Cover via LP Rounding")
    universe = {1, 2, 3, 4, 5}
    sets = {
        0: {1, 2, 3},
        1: {3, 4, 5},
        2: {1, 4},
        3: {2, 5}
    }
    costs = {0: 3, 1: 3, 2: 2, 3: 2}
    
    x, opt = set_cover_lp(universe, sets, costs)
    print(f"  LP solution: {x}")
    print(f"  LP optimal: {opt}")
    
    picked, cost = set_cover_lp_rounding(universe, sets, costs)
    print(f"  Rounded (f=2): {picked}, cost={cost}")
    
    picked2, cost2 = set_cover_primal_dual(universe, sets, costs)
    print(f"  Primal-Dual: {picked2}, cost={cost2}")
    
    # Vertex Cover via LP
    print("\n2. Vertex Cover via LP Rounding (2-approx)")
    graph = {
        0: {1: 1, 2: 1},
        1: {0: 1, 2: 1, 3: 1},
        2: {0: 1, 1: 1, 3: 1},
        3: {1: 1, 2: 1}
    }
    cover, cost = vertex_cover_lp_rounding(graph)
    print(f"  Graph: C4 (cycle of 4)")
    print(f"  LP-rounded cover: {cover}, size={cost}")
    print(f"  Optimal: {2}")


if __name__ == "__main__":
    demo_lp_algorithms()
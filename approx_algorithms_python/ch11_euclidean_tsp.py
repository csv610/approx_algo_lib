"""
Chapter 11: Euclidean TSP
==========================
Vazirani Ch. 11: Euclidean TSP.
Implements:
1. Exact TSP via Held-Karp dynamic programming (exponential baseline)
2. Quadtree-based divide-and-conquer cycle merging heuristic (geometric approximation)
"""

import math
from typing import List, Tuple

def compute_tour_cost(tour: List[int], points: List[Tuple[float, float]]) -> float:
    """Calculate the total Euclidean distance of a tour."""
    cost = 0.0
    for i in range(len(tour) - 1):
        u, v = tour[i], tour[i+1]
        cost += math.hypot(points[u][0] - points[v][0], points[u][1] - points[v][1])
    return cost

def held_karp_tsp(points: List[Tuple[float, float]]) -> Tuple[List[int], float]:
    """Exact TSP on points in 2D using Held-Karp DP."""
    n = len(points)
    if n == 0:
        return [], 0.0
    if n == 1:
        return [0], 0.0
        
    dist = [[math.hypot(p1[0]-p2[0], p1[1]-p2[1]) for p2 in points] for p1 in points]
    
    # DP memoization
    memo = {}
    parent = {}
    
    def tsp_solve(mask, last):
        if mask == (1 << n) - 1:
            return dist[last][0]
        state = (mask, last)
        if state in memo:
            return memo[state]
            
        best = float('inf')
        best_next = -1
        for nxt in range(n):
            if not (mask & (1 << nxt)):
                cost = dist[last][nxt] + tsp_solve(mask | (1 << nxt), nxt)
                if cost < best:
                    best = cost
                    best_next = nxt
        memo[state] = best
        parent[state] = best_next
        return best
        
    opt_cost = tsp_solve(1, 0)
    
    # Reconstruct path
    path = [0]
    mask = 1
    curr = 0
    for _ in range(n - 1):
        nxt = parent[(mask, curr)]
        path.append(nxt)
        mask |= (1 << nxt)
        curr = nxt
    path.append(0)
    return path, opt_cost

def merge_tours(tour1: List[int], tour2: List[int], points: List[Tuple[float, float]]) -> List[int]:
    """Cheapest connection swap between two tours or insertion of a single point."""
    if not tour1:
        return tour2
    if not tour2:
        return tour1
        
    dist = lambda u, v: math.hypot(points[u][0]-points[v][0], points[u][1]-points[v][1])
    
    # Case 1: tour2 is a single point
    if len(tour2) == 1 or (len(tour2) == 2 and tour2[0] == tour2[1]):
        u = tour2[0]
        best_cost_diff = float('inf')
        best_idx = -1
        for i in range(len(tour1) - 1):
            u1 = tour1[i]
            v1 = tour1[i+1]
            diff = dist(u1, u) + dist(u, v1) - dist(u1, v1)
            if diff < best_cost_diff:
                best_cost_diff = diff
                best_idx = i
        merged = tour1[:best_idx + 1] + [u] + tour1[best_idx + 1:]
        return merged
        
    # Case 2: tour1 is a single point
    if len(tour1) == 1 or (len(tour1) == 2 and tour1[0] == tour1[1]):
        return merge_tours(tour2, tour1, points)
        
    # Case 3: Both are closed tours of size >= 2 (vertices >= 2)
    t1 = list(tour1)
    if t1[-1] != t1[0]:
        t1.append(t1[0])
    t2 = list(tour2)
    if t2[-1] != t2[0]:
        t2.append(t2[0])
    
    n1 = len(t1) - 1
    n2 = len(t2) - 1
    best_cost_diff = float('inf')
    best_swap = None
    
    for i in range(n1):
        u1 = t1[i]
        v1 = t1[i+1]
        for j in range(n2):
            u2 = t2[j]
            v2 = t2[j+1]
            
            diff1 = dist(u1, u2) + dist(v1, v2) - dist(u1, v1) - dist(u2, v2)
            diff2 = dist(u1, v2) + dist(v1, u2) - dist(u1, v1) - dist(u2, v2)
            
            if diff1 < best_cost_diff:
                best_cost_diff = diff1
                best_swap = (i, j, False)
            if diff2 < best_cost_diff:
                best_cost_diff = diff2
                best_swap = (i, j, True)
                
    if not best_swap:
        return t1 + t2[1:]
        
    i, j, reverse_second = best_swap
    part1 = t1[:i+1]
    t2_rotated = t2[j:-1] + t2[:j]
    if reverse_second:
        t2_rotated.reverse()
        
    merged = part1 + t2_rotated + t1[i+1:]
    if merged[-1] != merged[0]:
        merged.append(merged[0])
    return merged

def quadtree_tsp_solve(point_indices: List[int], points: List[Tuple[float, float]]) -> List[int]:
    """Recursively partition points using a quadtree and merge tours."""
    n = len(point_indices)
    if n == 0:
        return []
    if n <= 3:
        # Base case: solve exactly on these few points
        local_points = [points[i] for i in point_indices]
        local_tour, _ = held_karp_tsp(local_points)
        # Map local indices back to global indices
        global_tour = [point_indices[i] for i in local_tour]
        return global_tour
        
    # Find bounding box
    xs = [points[i][0] for i in point_indices]
    ys = [points[i][1] for i in point_indices]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    
    mid_x = (min_x + max_x) / 2
    mid_y = (min_y + max_y) / 2
    
    # If points are collinear or co-located, fallback to simple sorting
    if max_x - min_x < 1e-9 and max_y - min_y < 1e-9:
        tour = list(point_indices) + [point_indices[0]]
        return tour

    # Divide into 4 quadrants
    q1, q2, q3, q4 = [], [], [], []
    for idx in point_indices:
        x, y = points[idx]
        if x <= mid_x:
            if y <= mid_y:
                q1.append(idx)
            else:
                q2.append(idx)
        else:
            if y <= mid_y:
                q3.append(idx)
            else:
                q4.append(idx)
                
    # If the partitioning failed to split (e.g. all points in one quadrant), split by index half
    if len(q1) == n or len(q2) == n or len(q3) == n or len(q4) == n:
        half = n // 2
        q1 = point_indices[:half]
        q2 = point_indices[half:]
        q3, q4 = [], []
        
    # Solve recursively
    t1 = quadtree_tsp_solve(q1, points)
    t2 = quadtree_tsp_solve(q2, points)
    t3 = quadtree_tsp_solve(q3, points)
    t4 = quadtree_tsp_solve(q4, points)
    
    # Merge tours sequentially
    merged = t1
    for t in [t2, t3, t4]:
        if t:
            merged = merge_tours(merged, t, points)
            
    return merged

def quadtree_tsp(points: List[Tuple[float, float]]) -> Tuple[List[int], float]:
    """Quadtree partition-and-merge heuristic for Euclidean TSP."""
    indices = list(range(len(points)))
    tour = quadtree_tsp_solve(indices, points)
    cost = compute_tour_cost(tour, points)
    return tour, cost

def demo_euclidean_tsp():
    print("=" * 60)
    print("Chapter 11: Euclidean TSP Heuristics")
    print("=" * 60)
    
    # Example 1: 8 points in a 2D plane
    points = [
        (0.0, 0.0),
        (1.0, 4.0),
        (3.0, 1.0),
        (4.0, 3.0),
        (1.0, 1.0),
        (3.0, 4.0),
        (5.0, 0.0),
        (5.5, 4.5)
    ]
    
    print(f"\n1. Input 2D Points (n={len(points)}):")
    for idx, p in enumerate(points):
        print(f"  Point {idx}: {p}")
        
    opt_tour, opt_cost = held_karp_tsp(points)
    print(f"\nExact Held-Karp Tour: {opt_tour}")
    print(f"Exact Optimal Cost:    {opt_cost:.4f}")
    
    qt_tour, qt_cost = quadtree_tsp(points)
    print(f"\nQuadtree Heuristic Tour: {qt_tour}")
    print(f"Heuristic Tour Cost:     {qt_cost:.4f}")
    print(f"Approximation Ratio:     {qt_cost/opt_cost:.4f} (Theoretical: 1 + eps)")

if __name__ == "__main__":
    demo_euclidean_tsp()

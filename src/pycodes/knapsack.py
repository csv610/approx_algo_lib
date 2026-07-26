"""
Chapter 8: Knapsack and FPTAS
===============================
Vazirani Ch. 8: Knapsack FPTAS via dynamic programming with rounding.
"""

from typing import List, Tuple, Dict
import math


def knapsack_dp(weights: List[int], values: List[int], capacity: int) -> Tuple[List[int], int]:
    """Exact 0/1 Knapsack via DP (pseudo-polynomial)."""
    n = len(weights)
    dp = [[0] * (capacity + 1) for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        for w in range(capacity + 1):
            if weights[i-1] <= w:
                dp[i][w] = max(dp[i-1][w], dp[i-1][w-weights[i-1]] + values[i-1])
            else:
                dp[i][w] = dp[i-1][w]
    
    # Reconstruct
    selected = []
    w = capacity
    for i in range(n, 0, -1):
        if dp[i][w] != dp[i-1][w]:
            selected.append(i-1)
            w -= weights[i-1]
    
    return selected, dp[n][capacity]


def knapsack_fptas(weights: List[int], values: List[int], capacity: int, 
                   epsilon: float) -> Tuple[List[int], int]:
    """
    FPTAS for 0/1 Knapsack (Vazirani, Section 8.2 / standard).
    
    Scale values: v_i' = floor(v_i * n / (epsilon * V_max))
    Then run DP on scaled values.
    
    Approximation: (1 - epsilon) * OPT
    Time: O(n^3 / epsilon)
    """
    n = len(weights)
    if n == 0:
        return [], 0
    
    V_max = max(values)
    if V_max == 0:
        return [], 0
    
    # Scaling factor
    K = (epsilon * V_max) / n
    if K == 0:
        K = 1
    
    # Scaled values
    scaled_values = [int(v / K) for v in values]
    
    # DP on scaled values
    # dp[v] = min weight to achieve value v
    max_scaled_val = sum(scaled_values)
    dp = [float('inf')] * (max_scaled_val + 1)
    dp[0] = 0
    
    # Track items
    choice = [[False] * (max_scaled_val + 1) for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        v = scaled_values[i-1]
        w = weights[i-1]
        for val in range(max_scaled_val, v - 1, -1):
            if dp[val - v] + w < dp[val]:
                dp[val] = dp[val - v] + w
                choice[i][val] = True
    
    # Find best feasible value
    best_val = 0
    for val in range(max_scaled_val + 1):
        if dp[val] <= capacity:
            best_val = val
    
    # Reconstruct
    selected = []
    val = best_val
    for i in range(n, 0, -1):
        if choice[i][val]:
            selected.append(i-1)
            val -= scaled_values[i-1]
    
    actual_value = sum(values[i] for i in selected)
    actual_weight = sum(weights[i] for i in selected)
    
    return selected, actual_value


def knapsack_greedy_ratio(weights: List[int], values: List[int], capacity: int) -> Tuple[List[int], int]:
    """Greedy by value/weight ratio (2-approx for fractional, bad for 0/1)."""
    n = len(weights)
    items = [(i, values[i] / weights[i] if weights[i] > 0 else float('inf')) 
             for i in range(n)]
    items.sort(key=lambda x: x[1], reverse=True)
    
    selected = []
    total_weight = 0
    total_value = 0
    
    for i, _ in items:
        if total_weight + weights[i] <= capacity:
            selected.append(i)
            total_weight += weights[i]
            total_value += values[i]
    
    return selected, total_value


def knapsack_unbounded_dp(weights: List[int], values: List[int], capacity: int) -> Tuple[List[int], int]:
    """Unbounded Knapsack (unlimited copies) - exact DP."""
    n = len(weights)
    dp = [0] * (capacity + 1)
    choice = [-1] * (capacity + 1)
    
    for w in range(1, capacity + 1):
        for i in range(n):
            if weights[i] <= w and dp[w - weights[i]] + values[i] > dp[w]:
                dp[w] = dp[w - weights[i]] + values[i]
                choice[w] = i
    
    # Reconstruct
    selected = []
    w = capacity
    while w > 0 and choice[w] != -1:
        selected.append(choice[w])
        w -= weights[choice[w]]
    
    return selected, dp[capacity]


def demo_knapsack():
    print("=" * 60)
    print("Chapter 8: Knapsack FPTAS")
    print("=" * 60)
    
    # Example from textbook
    weights = [10, 20, 30, 40, 50]
    values = [60, 100, 120, 200, 250]
    capacity = 100
    
    print(f"\nInstance: n={len(weights)}, capacity={capacity}")
    print(f"Weights: {weights}")
    print(f"Values:  {values}")
    
    # Exact DP
    sel_exact, val_exact = knapsack_dp(weights, values, capacity)
    print(f"\nExact DP: items={sel_exact}, value={val_exact}")
    
    # FPTAS for different epsilon
    for eps in [0.5, 0.25, 0.1, 0.05, 0.01]:
        sel, val = knapsack_fptas(weights, values, capacity, eps)
        ratio = val / val_exact if val_exact > 0 else 0
        print(f"  FPTAS eps={eps:.2f}: value={val}, ratio={ratio:.4f} (bound={1-eps:.4f})")
    
    # Another example
    print("\n--- Another Example ---")
    weights2 = [2, 3, 4, 5]
    values2 = [3, 4, 5, 6]
    capacity2 = 8
    
    sel_exact, val_exact = knapsack_dp(weights2, values2, capacity2)
    print(f"Exact: {sel_exact}, value={val_exact}")
    
    for eps in [0.2, 0.1]:
        sel, val = knapsack_fptas(weights2, values2, capacity2, eps)
        print(f"  FPTAS eps={eps}: {sel}, value={val}, ratio={val/val_exact:.4f}")
    
    # Unbounded knapsack
    print("\n--- Unbounded Knapsack ---")
    sel_ub, val_ub = knapsack_unbounded_dp(weights, values, capacity)
    print(f"Unbounded exact: items={sel_ub}, value={val_ub}")


if __name__ == "__main__":
    demo_knapsack()

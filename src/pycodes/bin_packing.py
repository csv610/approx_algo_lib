"""
Chapter 9: Bin Packing
======================
Vazirani Ch. 9: Bin Packing approximation algorithms.
Implements:
1. Next-Fit (2-approx online)
2. First-Fit (1.7-approx online)
3. First-Fit Decreasing (11/9-approx offline)
4. Asymptotic PTAS (APTAS) via linear grouping and DP for large items.
"""

import math
from typing import List, Tuple

def next_fit(items: List[float], capacity: float = 1.0) -> List[List[float]]:
    """Next-Fit (NF) 2-approximation online bin packing algorithm."""
    bins = []
    current_bin = []
    current_weight = 0.0
    
    for item in items:
        if current_weight + item <= capacity:
            current_bin.append(item)
            current_weight += item
        else:
            if current_bin:
                bins.append(current_bin)
            current_bin = [item]
            current_weight = item
            
    if current_bin:
        bins.append(current_bin)
    return bins

def first_fit(items: List[float], capacity: float = 1.0) -> List[List[float]]:
    """First-Fit (FF) 1.7-approximation online bin packing algorithm."""
    bins = []
    
    for item in items:
        placed = False
        for b in bins:
            if sum(b) + item <= capacity:
                b.append(item)
                placed = True
                break
        if not placed:
            bins.append([item])
            
    return bins

def first_fit_decreasing(items: List[float], capacity: float = 1.0) -> List[List[float]]:
    """First-Fit Decreasing (FFD) 11/9-approximation offline bin packing algorithm."""
    # Sort items in descending order
    sorted_items = sorted(items, reverse=True)
    return first_fit(sorted_items, capacity)

def generate_configurations(sizes: List[float], cap: float = 1.0) -> List[Tuple[int, ...]]:
    """Generate all valid item count configurations that fit in a bin."""
    configs = []
    n = len(sizes)
    
    def backtrack(idx, current_conf, remaining_cap):
        if idx == n:
            if sum(current_conf) > 0:
                configs.append(tuple(current_conf))
            return
        max_count = int(remaining_cap // sizes[idx])
        for count in range(max_count + 1):
            current_conf.append(count)
            backtrack(idx + 1, current_conf, remaining_cap - count * sizes[idx])
            current_conf.pop()
            
    backtrack(0, [], cap)
    return configs

def pack_large_dp(counts: Tuple[int, ...], configs: List[Tuple[int, ...]], sizes: List[float]) -> List[List[float]]:
    """Solve the rounded large items bin packing instance exactly using DP with memoization."""
    memo = {}
    
    def solve(state: Tuple[int, ...]) -> Tuple[int, List[List[float]]]:
        if sum(state) == 0:
            return 0, []
        if state in memo:
            return memo[state]
            
        best_val = float('inf')
        best_bins = []
        
        for conf in configs:
            if all(state[i] >= conf[i] for i in range(len(state))):
                next_state = tuple(state[i] - conf[i] for i in range(len(state)))
                val, bins = solve(next_state)
                if 1 + val < best_val:
                    best_val = 1 + val
                    # Create a new bin with items from conf
                    new_bin = []
                    for i, count in enumerate(conf):
                        new_bin.extend([sizes[i]] * count)
                    # We store the bins as list of lists
                    best_bins = [new_bin] + [list(b) for b in bins]
                    
        memo[state] = (best_val, best_bins)
        return memo[state]
        
    return solve(counts)[1]

def bin_packing_aptas(items: List[float], eps: float = 0.3, capacity: float = 1.0) -> List[List[float]]:
    """
    Asymptotic PTAS for Bin Packing via linear grouping.
    Packs items into at most (1 + 2*eps)*OPT + O(1/eps^2) bins.
    """
    # 1. Separate large and small items
    large_items = [x for x in items if x >= eps]
    small_items = [x for x in items if x < eps]
    
    if not large_items:
        # If no large items, just pack small items using First-Fit
        return first_fit(small_items, capacity)
        
    # 2. Linear Grouping on large items
    large_items.sort()
    n_large = len(large_items)
    k = int(1.0 / (eps ** 2))
    q = n_large // k  # size of each group
    
    # Map from original large item index to rounded size
    rounded_large = []
    
    if q == 0 or k == 0:
        # If there are too few items to group, solve the exact large instance
        rounded_large = list(large_items)
    else:
        # Group items and round up to the maximum size of the group
        groups = []
        for i in range(k):
            start = i * q
            end = (i + 1) * q if i < k - 1 else n_large
            group = large_items[start:end]
            groups.append(group)
            
        # Round up size to the max of the group
        for i, group in enumerate(groups):
            max_size = max(group)
            rounded_large.extend([max_size] * len(group))
            
    # 3. Solve the rounded instance of large items via DP
    # Identify distinct sizes and their counts
    distinct_sizes = sorted(list(set(rounded_large)))
    counts = tuple(rounded_large.count(s) for s in distinct_sizes)
    
    # Generate valid configurations
    configs = generate_configurations(distinct_sizes, capacity)
    
    # Solve DP
    large_bins = pack_large_dp(counts, configs, distinct_sizes)
    
    # Map rounded items back to original large items (in sorted order)
    # Since any packing of rounded items is valid for smaller original items:
    # We sort the items inside large_bins, and replace them with the sorted original large_items.
    # To do this safely: we flatten the bins, replace values, and reshape back.
    flattened_bins_items = []
    for b in large_bins:
        flattened_bins_items.extend(b)
    flattened_bins_items.sort()
    
    # Map index
    item_map = {}
    for r_val, orig_val in zip(flattened_bins_items, large_items):
        item_map[r_val] = item_map.get(r_val, []) + [orig_val]
        
    final_large_bins = []
    for b in large_bins:
        new_bin = []
        for item in b:
            # Pop one original value that corresponds to this rounded size
            orig_val = item_map[item].pop(0)
            new_bin.append(orig_val)
        final_large_bins.append(new_bin)
        
    # 4. Pack small items using First-Fit into the existing bins or open new ones
    for item in small_items:
        placed = False
        for b in final_large_bins:
            if sum(b) + item <= capacity:
                b.append(item)
                placed = True
                break
        if not placed:
            final_large_bins.append([item])
            
    return final_large_bins

def demo_bin_packing():
    print("=" * 60)
    print("Chapter 9: Bin Packing Algorithms")
    print("=" * 60)
    
    # Example 1: Basic items
    items1 = [0.2, 0.5, 0.4, 0.7, 0.1, 0.3, 0.8]
    print(f"\n1. Input Items: {items1}")
    print(f"  Next-Fit (NF):           {next_fit(items1)} (bins: {len(next_fit(items1))})")
    print(f"  First-Fit (FF):          {first_fit(items1)} (bins: {len(first_fit(items1))})")
    print(f"  First-Fit Decreasing:    {first_fit_decreasing(items1)} (bins: {len(first_fit_decreasing(items1))})")
    print(f"  APTAS (eps=0.4):         {bin_packing_aptas(items1, eps=0.4)} (bins: {len(bin_packing_aptas(items1, eps=0.4))})")
    
    # Example 2: Larger random-like instance
    items2 = [0.15, 0.33, 0.45, 0.12, 0.61, 0.38, 0.49, 0.52, 0.23, 0.29, 0.41, 0.19, 0.31, 0.27, 0.55, 0.48, 0.35, 0.11]
    print(f"\n2. Larger Instance (n={len(items2)}):")
    print(f"  Next-Fit bins:           {len(next_fit(items2))}")
    print(f"  First-Fit bins:          {len(first_fit(items2))}")
    print(f"  First-Fit Decreasing:    {len(first_fit_decreasing(items2))}")
    print(f"  APTAS (eps=0.35) bins:   {len(bin_packing_aptas(items2, eps=0.35))}")

if __name__ == "__main__":
    demo_bin_packing()

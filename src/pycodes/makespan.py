"""
Chapter 10: Minimum Makespan Scheduling
=======================================
Vazirani Ch. 10: Minimum Makespan Scheduling.
Implements:
1. Graham's List Scheduling (online 2-approx)
2. Longest Processing Time (LPT) scheduling (offline 4/3-approx)
3. Hochbaum-Shmoys PTAS (Polynomial-Time Approximation Scheme)
"""

import math
from typing import List, Tuple
from bin_packing import generate_configurations, pack_large_dp

def list_scheduling(jobs: List[float], m: int) -> List[List[float]]:
    """Graham's List Scheduling algorithm (online 2-approximation)."""
    schedule = [[] for _ in range(m)]
    for job in jobs:
        machine_idx = min(range(m), key=lambda i: sum(schedule[i]))
        schedule[machine_idx].append(job)
    return schedule

def lpt_scheduling(jobs: List[float], m: int) -> List[List[float]]:
    """Longest Processing Time (LPT) scheduling algorithm (offline 4/3-approximation)."""
    sorted_jobs = sorted(jobs, reverse=True)
    return list_scheduling(sorted_jobs, m)

def check_schedule(jobs: List[float], m: int, T: float, eps: float) -> Tuple[bool, List[List[float]]]:
    """Test if the jobs can be scheduled on m machines within makespan (1+eps)*T using dual approximation."""
    large_jobs = [p for p in jobs if p > eps * T]
    small_jobs = [p for p in jobs if p <= eps * T]
    
    if not large_jobs:
        schedule = [[] for _ in range(m)]
        for job in small_jobs:
            machine_idx = min(range(m), key=lambda i: sum(schedule[i]))
            schedule[machine_idx].append(job)
        max_load = max(sum(sch) for sch in schedule)
        return max_load <= (1 + eps) * T, schedule

    # Round down large jobs
    delta = (eps ** 2) * T
    rounded_large = []
    orig_large_sorted = sorted(large_jobs)
    for p in orig_large_sorted:
        val = math.floor(p / delta) * delta
        val = max(val, eps * T)
        rounded_large.append(val)
        
    distinct_sizes = sorted(list(set(rounded_large)))
    counts = tuple(rounded_large.count(s) for s in distinct_sizes)
    configs = generate_configurations(distinct_sizes, T)
    
    large_bins = pack_large_dp(counts, configs, distinct_sizes)
    
    if len(large_bins) > m:
        return False, [[] for _ in range(m)]
        
    # Map rounded large jobs back to original large jobs
    flattened_bins_items = []
    for b in large_bins:
        flattened_bins_items.extend(b)
    flattened_bins_items.sort()
    
    item_map = {s: [] for s in distinct_sizes}
    for r_val, orig_val in zip(flattened_bins_items, orig_large_sorted):
        item_map[r_val].append(orig_val)
        
    schedule = [[] for _ in range(m)]
    for idx, b in enumerate(large_bins):
        for item in b:
            orig_val = item_map[item].pop(0)
            schedule[idx].append(orig_val)
            
    # Place small jobs greedily
    for job in small_jobs:
        machine_idx = min(range(m), key=lambda i: sum(schedule[i]))
        schedule[machine_idx].append(job)
        
    max_load = max(sum(sch) for sch in schedule)
    return max_load <= (1 + eps) * T, schedule

def makespan_ptas(jobs: List[float], m: int, eps: float = 0.25) -> List[List[float]]:
    """Hochbaum-Shmoys PTAS for Minimum Makespan Scheduling."""
    lb = max(max(jobs), sum(jobs) / m)
    ub = max(max(jobs), 2 * sum(jobs) / m)
    
    best_schedule = None
    for _ in range(15):
        mid = (lb + ub) / 2
        success, schedule = check_schedule(jobs, m, mid, eps)
        if success:
            best_schedule = schedule
            ub = mid
        else:
            lb = mid
            
    if not best_schedule:
        _, best_schedule = check_schedule(jobs, m, ub, eps)
    return best_schedule

def demo_makespan():
    print("=" * 60)
    print("Chapter 10: Minimum Makespan Scheduling")
    print("=" * 60)
    
    # 5 jobs, 2 machines
    jobs1 = [2.0, 3.0, 4.0, 6.0, 2.0]
    m1 = 2
    print(f"\n1. Jobs: {jobs1} on {m1} machines")
    
    sched_list = list_scheduling(jobs1, m1)
    print(f"  List Scheduling:      {sched_list} (makespan: {max(sum(s) for s in sched_list)})")
    
    sched_lpt = lpt_scheduling(jobs1, m1)
    print(f"  LPT Heuristic:        {sched_lpt} (makespan: {max(sum(s) for s in sched_lpt)})")
    
    sched_ptas = makespan_ptas(jobs1, m1, eps=0.25)
    print(f"  PTAS (eps=0.25):      {sched_ptas} (makespan: {max(sum(s) for s in sched_ptas)})")
    
    # Larger job list
    jobs2 = [1.2, 2.5, 3.1, 4.0, 1.8, 2.2, 5.0, 3.5, 0.9, 1.6]
    m2 = 3
    print(f"\n2. Larger Instance (n={len(jobs2)}) on {m2} machines:")
    
    sched_list2 = list_scheduling(jobs2, m2)
    print(f"  List Scheduling makespan:   {max(sum(s) for s in sched_list2)}")
    
    sched_lpt2 = lpt_scheduling(jobs2, m2)
    print(f"  LPT Heuristic makespan:     {max(sum(s) for s in sched_lpt2)}")
    
    sched_ptas2 = makespan_ptas(jobs2, m2, eps=0.2)
    print(f"  PTAS (eps=0.2) makespan:    {max(sum(s) for s in sched_ptas2)}")

if __name__ == "__main__":
    demo_makespan()

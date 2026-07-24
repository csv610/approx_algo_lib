"""
Chapter 7: Shortest Superstring
===============================
Vazirani Ch. 7: Shortest Superstring approximation algorithms.
Implements the 4-approximation (cycle cover concatenation) and 
the 3-approximation (modified greedy cycle-cover merge).
"""

from typing import List, Tuple, Set

def compute_overlap(s1: str, s2: str) -> int:
    """Length of the longest suffix of s1 that is a prefix of s2."""
    max_len = min(len(s1), len(s2))
    for l in range(max_len, 0, -1):
        if s1.endswith(s2[:l]):
            return l
    return 0

def preprocess_substrings(strings: List[str]) -> List[str]:
    """Filter out strings that are substrings of other strings."""
    sorted_strings = sorted(strings, key=len, reverse=True)
    filtered = []
    for s in sorted_strings:
        if not any(s in other for other in filtered):
            filtered.append(s)
    return filtered

def find_minimum_cycle_cover(cost_matrix: List[List[float]]) -> Tuple[List[int], float]:
    """Find a minimum-weight cycle cover using simple backtracking search (suitable for small n)."""
    n = len(cost_matrix)
    best_cost = float('inf')
    best_perm = []
    
    def backtrack(curr, perm, visited, cost):
        nonlocal best_cost, best_perm
        if cost >= best_cost:
            return
        if curr == n:
            best_cost = cost
            best_perm = list(perm)
            return
        for next_val in range(n):
            if not visited[next_val]:
                visited[next_val] = True
                perm.append(next_val)
                backtrack(curr + 1, perm, visited, cost + cost_matrix[curr][next_val])
                perm.pop()
                visited[next_val] = False

    backtrack(0, [], [False] * n, 0.0)
    return best_perm, best_cost

def extract_cycles(perm: List[int]) -> List[List[int]]:
    """Extract cycles from the permutation list."""
    n = len(perm)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if not visited[i]:
            cycle = []
            curr = i
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = perm[curr]
            cycles.append(cycle)
    return cycles

def greedy_superstring(strings: List[str]) -> str:
    """Greedy algorithm for Shortest Superstring."""
    T = list(strings)
    while len(T) > 1:
        max_ov = -1
        best_i, best_j = -1, -1
        for i in range(len(T)):
            for j in range(len(T)):
                if i != j:
                    ov = compute_overlap(T[i], T[j])
                    if ov > max_ov:
                        max_ov = ov
                        best_i, best_j = i, j
        
        # Merge best_i and best_j
        s_i = T[best_i]
        s_j = T[best_j]
        merged = s_i + s_j[max_ov:]
        
        # Remove and append
        if best_i < best_j:
            T.pop(best_j)
            T.pop(best_i)
        else:
            T.pop(best_i)
            T.pop(best_j)
        T.append(merged)
        
    return T[0] if T else ""

def cycle_cover_to_strings(cycles: List[List[int]], strings: List[str]) -> List[str]:
    """Merge the strings in each cycle to form one representative string per cycle."""
    cycle_strings = []
    for cycle in cycles:
        # If the cycle is a single node, its representative is the string itself
        if len(cycle) == 1:
            cycle_strings.append(strings[cycle[0]])
            continue
        # Otherwise, merge the strings sequentially along the cycle
        merged = strings[cycle[0]]
        for idx in range(len(cycle) - 1):
            u, v = cycle[idx], cycle[idx + 1]
            ov = compute_overlap(strings[u], strings[v])
            merged += strings[v][ov:]
        # Wrap around: overlap between last and first
        last, first = cycle[-1], cycle[0]
        ov = compute_overlap(strings[last], strings[first])
        # Note: we do not double-include the first string, but the cycle cover math
        # accounts for this overlap.
        cycle_strings.append(merged)
    return cycle_strings

def shortest_superstring_4approx(strings: List[str]) -> str:
    """4-approximation algorithm for Shortest Superstring (Blum et al.)."""
    cleaned = preprocess_substrings(strings)
    if not cleaned:
        return ""
    if len(cleaned) == 1:
        return cleaned[0]
        
    n = len(cleaned)
    # Build cost matrix: cost(u, v) = |v| - overlap(u, v)
    cost_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cost_matrix[i][j] = len(cleaned[j]) - compute_overlap(cleaned[i], cleaned[j])
            
    perm, _ = find_minimum_cycle_cover(cost_matrix)
    cycles = extract_cycles(perm)
    cycle_strings = cycle_cover_to_strings(cycles, cleaned)
    
    # Concat all cycle representative strings
    return "".join(cycle_strings)

def shortest_superstring_3approx(strings: List[str]) -> str:
    """3-approximation algorithm for Shortest Superstring (Modified Greedy)."""
    cleaned = preprocess_substrings(strings)
    if not cleaned:
        return ""
    if len(cleaned) == 1:
        return cleaned[0]
        
    n = len(cleaned)
    cost_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cost_matrix[i][j] = len(cleaned[j]) - compute_overlap(cleaned[i], cleaned[j])
            
    perm, _ = find_minimum_cycle_cover(cost_matrix)
    cycles = extract_cycles(perm)
    cycle_strings = cycle_cover_to_strings(cycles, cleaned)
    
    # Run greedy algorithm on the cycle representative strings
    return greedy_superstring(cycle_strings)

def demo_shortest_superstring():
    print("=" * 60)
    print("Chapter 7: Shortest Common Superstring")
    print("=" * 60)
    
    # Test case 1
    S1 = ["abc", "bcd", "cde", "def"]
    print(f"\n1. Input Strings: {S1}")
    print(f"  Greedy Superstring:  {greedy_superstring(S1)}")
    print(f"  4-Approx Superstring: {shortest_superstring_4approx(S1)}")
    print(f"  3-Approx Superstring: {shortest_superstring_3approx(S1)}")
    
    # Test case 2: DNA sequencing fragments
    S2 = ["CATG", "ATGT", "TGTA", "GTAC", "TACA"]
    print(f"\n2. DNA Fragments: {S2}")
    print(f"  Greedy Superstring:  {greedy_superstring(S2)}")
    print(f"  4-Approx Superstring: {shortest_superstring_4approx(S2)}")
    print(f"  3-Approx Superstring: {shortest_superstring_3approx(S2)}")

if __name__ == "__main__":
    demo_shortest_superstring()

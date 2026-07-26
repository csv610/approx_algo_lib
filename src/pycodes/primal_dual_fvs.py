"""
Chapter 23: Feedback Vertex Set via Primal-Dual
================================================
Vazirani Ch. 23 / Undirected FVS primal-dual 2-approximation.
Implements:
1. Bafna-Berman-Fujito primal-dual algorithm for undirected FVS.
2. Degree-based weight reduction over active cycles.
3. Reverse-delete pruning to guarantee a factor 2 approximation.
"""

from typing import List, Tuple, Dict, Set

def find_any_cycle(vertices: Set[int], edges: List[Tuple[int, int]]) -> List[int]:
    """Finds a simple cycle in the graph if one exists, using recursive DFS."""
    adj = {v: [] for v in vertices}
    for u, v in edges:
        if u in adj and v in adj:
            adj[u].append(v)
            adj[v].append(u)
            
    visited = set()
    parent = {}
    
    def dfs(node: int, p: int) -> List[int]:
        visited.add(node)
        parent[node] = p
        for neighbor in adj[node]:
            if neighbor != p:
                if neighbor in visited:
                    # Cycle found: trace back from node to neighbor
                    cycle = []
                    curr = node
                    while curr != neighbor:
                        cycle.append(curr)
                        curr = parent[curr]
                    cycle.append(neighbor)
                    return cycle
                else:
                    res = dfs(neighbor, node)
                    if res:
                        return res
        return None
        
    for start in vertices:
        if start not in visited:
            res = dfs(start, -1)
            if res:
                return res
    return []

def is_acyclic(vertices: Set[int], edges: List[Tuple[int, int]], removed: Set[int]) -> bool:
    """Returns True if the graph after removing vertices is acyclic (a forest)."""
    remaining_v = vertices - removed
    remaining_e = [(u, v) for u, v in edges if u in remaining_v and v in remaining_v]
    return len(find_any_cycle(remaining_v, remaining_e)) == 0

def primal_dual_fvs(
    vertices: List[int],
    edges: List[Tuple[int, int]],
    weights: Dict[int, float]
) -> List[int]:
    """
    Bafna-Berman-Fujito primal-dual 2-approximation for undirected FVS.
    """
    w = dict(weights)
    S = []
    
    active_vertices = set(vertices)
    active_edges = list(edges)
    
    while True:
        # 1. Clean the graph: recursively remove vertices of degree <= 1
        while True:
            deg = {v: 0 for v in active_vertices}
            for u, v in active_edges:
                deg[u] += 1
                deg[v] += 1
            to_remove = [v for v in active_vertices if deg[v] <= 1]
            if not to_remove:
                break
            for v in to_remove:
                active_vertices.remove(v)
                active_edges = [(u, x) for u, x in active_edges if u != v and x != v]
                
        if not active_vertices:
            break
            
        # 2. Find a cycle in the active graph
        cycle = find_any_cycle(active_vertices, active_edges)
        if not cycle:
            break
            
        # Compute current degrees in the active graph
        deg = {v: 0 for v in active_vertices}
        for u, v in active_edges:
            deg[u] += 1
            deg[v] += 1
            
        # 3. Grow dual variable on this cycle C
        # Rate of slack decrease for v in C is d(v) - 1
        delta = float('inf')
        for v in cycle:
            rate = deg[v] - 1
            if rate > 0:
                val = w[v] / rate
                if val < delta:
                    delta = val
                    
        # Update weights and identify the tight vertex
        tight_vertex = None
        for v in cycle:
            rate = deg[v] - 1
            w[v] -= delta * rate
            if w[v] <= 1e-9 and tight_vertex is None:
                tight_vertex = v
                
        if tight_vertex is None:
            tight_vertex = cycle[0]
            
        S.append(tight_vertex)
        active_vertices.remove(tight_vertex)
        active_edges = [(u, v) for u, v in active_edges if u != tight_vertex and v != tight_vertex]
        
    # 4. Reverse-delete pruning phase
    pruned = list(S)
    for v in reversed(S):
        pruned.remove(v)
        if not is_acyclic(set(vertices), edges, set(pruned)):
            pruned.append(v)
            
    return pruned

def demo_primal_dual_fvs():
    print("=" * 60)
    print("Chapter 23: Feedback Vertex Set via Primal-Dual")
    print("=" * 60)
    
    # Example 1: Intersecting Cycles (Shared Center Vertex)
    # Cycle 1: 0-1-2-0
    # Cycle 2: 2-3-4-2
    # Vertex 2 is shared. Optimal FVS is {2} (cost 3)
    vertices = [0, 1, 2, 3, 4]
    edges = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]
    weights = {0: 10.0, 1: 10.0, 2: 3.0, 3: 10.0, 4: 10.0}
    
    print("\n1. Intersecting Cycles (shared center node 2):")
    print(f"  Vertices: {vertices}")
    print(f"  Edges:    {edges}")
    print(f"  Weights:  {weights}")
    
    fvs = primal_dual_fvs(vertices, edges, weights)
    print(f"  Selected FVS: {fvs}")
    print(f"  Total Cost:   {sum(weights[v] for v in fvs):.2f}")
    
    # Example 2: Bipartite Graph K_3,3 (optimal is 4, approx within 2)
    # K_3,3: Left: 0,1,2, Right: 3,4,5
    # All weights 1
    v2 = list(range(6))
    e2 = [(u, v) for u in [0, 1, 2] for v in [3, 4, 5]]
    w2 = {i: 1.0 for i in v2}
    
    print("\n2. Bipartite Graph K_3,3 (unit weights):")
    fvs2 = primal_dual_fvs(v2, e2, w2)
    print(f"  Selected FVS: {fvs2}")
    print(f"  Total Cost:   {sum(w2[v] for v in fvs2):.2f}")

if __name__ == "__main__":
    demo_primal_dual_fvs()

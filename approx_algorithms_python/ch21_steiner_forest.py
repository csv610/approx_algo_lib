"""
Chapter 21: Steiner Forest
==========================
Vazirani Ch. 21: Steiner Forest via Primal-Dual.
Implements:
1. Agrawal-Klein-Ravi (AKR) primal-dual algorithm.
2. Reverse-delete pruning phase.
"""

from typing import List, Tuple, Dict, Set

def find_components(vertices: List[int], forest_edges: List[Tuple[int, int]]) -> List[Set[int]]:
    """Find connected components in G=(V, F)."""
    adj = {v: [] for v in vertices}
    for u, v in forest_edges:
        adj[u].append(v)
        adj[v].append(u)
        
    visited = set()
    components = []
    
    for v in vertices:
        if v not in visited:
            comp = set()
            stack = [v]
            while stack:
                curr = stack.pop()
                if curr not in comp:
                    comp.add(curr)
                    visited.add(curr)
                    for nxt in adj[curr]:
                        if nxt not in visited:
                            stack.append(nxt)
            components.append(comp)
    return components

def is_connected(
    vertices: List[int],
    forest_edges: List[Tuple[int, int]],
    s: int,
    t: int
) -> bool:
    """Check if s and t are connected in G=(V, F)."""
    adj = {v: [] for v in vertices}
    for u, v in forest_edges:
        adj[u].append(v)
        adj[v].append(u)
        
    visited = {s}
    stack = [s]
    while stack:
        curr = stack.pop()
        if curr == t:
            return True
        for nxt in adj[curr]:
            if nxt not in visited:
                visited.add(nxt)
                stack.append(nxt)
    return False

def steiner_forest_primal_dual(
    vertices: List[int],
    edges: List[Tuple[int, int, float]],
    pairs: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """
    AKR Primal-Dual 2-approximation algorithm for Steiner Forest.
    
    Args:
        vertices: List of vertex indices.
        edges: List of edges as (u, v, cost).
        pairs: List of terminal pairs (s_i, t_i) to connect.
    """
    # Active/tight edge selection phase
    chosen_edges = []
    L = {i: 0.0 for i in range(len(edges))} # accumulated dual sum on edges
    
    while True:
        # 1. Find connected components of G=(V, chosen_edges)
        components = find_components(vertices, chosen_edges)
        
        # Map node to component index
        node_to_comp = {}
        for idx, comp in enumerate(components):
            for node in comp:
                node_to_comp[node] = idx
                
        # 2. Identify active components
        # A component is active if it contains exactly one of the terminals of a pair
        # that are not yet connected.
        active_indices = set()
        for s, t in pairs:
            c_s = node_to_comp[s]
            c_t = node_to_comp[t]
            if c_s != c_t:
                active_indices.add(c_s)
                active_indices.add(c_t)
                
        # If no active components exist, we have connected all pairs
        if not active_indices:
            break
            
        # 3. Find next edge to become tight
        # Compute rates of dual increase for each edge
        best_delta = float('inf')
        best_edge_indices = []
        
        edge_rates = []
        for i, (u, v, cost) in enumerate(edges):
            if (u, v) in chosen_edges or (v, u) in chosen_edges:
                edge_rates.append(0)
                continue
                
            c_u = node_to_comp[u]
            c_v = node_to_comp[v]
            
            if c_u == c_v:
                # Inside same component, no dual contribution
                edge_rates.append(0)
                continue
                
            rate = 0
            if c_u in active_indices:
                rate += 1
            if c_v in active_indices:
                rate += 1
                
            edge_rates.append(rate)
            if rate > 0:
                delta = (cost - L[i]) / rate
                if delta < best_delta:
                    best_delta = delta
                    best_edge_indices = [i]
                elif abs(delta - best_delta) < 1e-9:
                    best_edge_indices.append(i)
                    
        # Update dual values on edges and add tight edges
        for i in range(len(edges)):
            if edge_rates[i] > 0:
                L[i] += best_delta * edge_rates[i]
                
        for idx in best_edge_indices:
            u, v, _ = edges[idx]
            chosen_edges.append((u, v))
            
    # 4. Reverse-delete pruning phase
    pruned_edges = list(chosen_edges)
    for edge in reversed(chosen_edges):
        pruned_edges.remove(edge)
        # Check if all pairs are still connected
        still_connected = True
        for s, t in pairs:
            if not is_connected(vertices, pruned_edges, s, t):
                still_connected = False
                break
        if not still_connected:
            # We need this edge, put it back
            pruned_edges.append(edge)
            
    return pruned_edges

def demo_steiner_forest():
    print("=" * 60)
    print("Chapter 21: Steiner Forest via Primal-Dual")
    print("=" * 60)
    
    # Example graph: path-like with cross-edges
    vertices = [0, 1, 2, 3, 4]
    edges = [
        (0, 1, 1.0),
        (1, 2, 2.0),
        (2, 3, 1.0),
        (3, 4, 3.0),
        (0, 4, 10.0) # detour edge
    ]
    
    # Connect 0 to 2, and 2 to 4
    pairs = [(0, 2), (2, 4)]
    
    print(f"\nGraph Vertices: {vertices}")
    print(f"Graph Edges:    {edges}")
    print(f"Terminal Pairs: {pairs}")
    
    forest = steiner_forest_primal_dual(vertices, edges, pairs)
    forest_cost = sum(cost for u, v, cost in edges if (u, v) in forest or (v, u) in forest)
    
    print(f"\nSelected Steiner Forest Edges: {forest}")
    print(f"Forest Total Cost:             {forest_cost:.2f}")

if __name__ == "__main__":
    demo_steiner_forest()

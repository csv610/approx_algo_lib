"""
Chapter 6: Feedback Vertex Set
===============================
Vazirani Ch. 6: Feedback Vertex Set 2-approximation in undirected weighted graphs.
Using local-ratio / layering technique.
"""

from typing import Dict, List, Set, Tuple

# Type aliases
Multigraph = Dict[int, List[int]]  # adjacency list allowing duplicate neighbors (multi-edges)

def connected_components_count(graph: Multigraph) -> int:
    """Count the number of connected components in a multigraph."""
    visited = set()
    count = 0
    for v in graph:
        if v not in visited:
            count += 1
            q = [v]
            visited.add(v)
            while q:
                curr = q.pop(0)
                for nxt in graph[curr]:
                    if nxt in graph and nxt not in visited:
                        visited.add(nxt)
                        q.append(nxt)
    return count

def cyclomatic_number(graph: Multigraph) -> int:
    """Calculate the cyclomatic number of the multigraph: m - n + cc."""
    n = len(graph)
    if n == 0:
        return 0
    # Sum degree and divide by 2 for edges
    m = sum(len(graph[v]) for v in graph) // 2
    cc = connected_components_count(graph)
    return m - n + cc

def is_feedback_vertex_set(graph: Multigraph, fvs: Set[int]) -> bool:
    """Check if the set fvs is a feedback vertex set of the graph (i.e. graph \\ fvs is acyclic)."""
    remaining = set(graph.keys()) - fvs
    visited = set()
    
    for start in remaining:
        if start not in visited:
            # BFS with parent tracking
            q = [(start, -1)]
            visited.add(start)
            while q:
                curr, parent = q.pop(0)
                neighbors = [x for x in graph[curr] if x in remaining]
                # If there are duplicate neighbors, it is a cycle of length 2
                if len(neighbors) != len(set(neighbors)):
                    return False
                for nxt in neighbors:
                    if nxt == parent:
                        continue
                    if nxt in visited:
                        return False  # Cycle detected
                    visited.add(nxt)
                    q.append((nxt, curr))
    return True

def feedback_vertex_set_approx(graph: Multigraph, weights: Dict[int, float]) -> Set[int]:
    """
    2-approximation for undirected weighted Feedback Vertex Set.
    Based on Vazirani's Chapter 6 layering/local ratio algorithm.
    """
    # Create a local copy of the graph
    g = {v: list(neighbors) for v, neighbors in graph.items()}
    
    # 0. If any vertex has weight <= 0, remove and recurse
    for v in list(g.keys()):
        if weights.get(v, 0.0) <= 1e-9:
            g_sub = {x: [y for y in neighbors if y != v] for x, neighbors in g.items() if x != v}
            sub_fvs = feedback_vertex_set_approx(g_sub, weights)
            sub_fvs.add(v)
            
            # Pruning before returning
            final_fvs = set(sub_fvs)
            for x in list(final_fvs):
                if is_feedback_vertex_set(graph, final_fvs - {x}):
                    final_fvs.remove(x)
            return final_fvs
    
    # 1. Clean the graph (remove degree <= 1 and contract degree 2)
    changed = True
    while changed:
        changed = False
        # Remove degree <= 1 vertices
        to_remove = [v for v in g if len(g[v]) <= 1]
        if to_remove:
            for v in to_remove:
                for nxt in g[v]:
                    if nxt in g:
                        while v in g[nxt]:
                            g[nxt].remove(v)
                del g[v]
            changed = True
            continue
            
        # Contract degree 2 vertices (only if neighbors are distinct)
        deg2_vertices = [v for v in g if len(g[v]) == 2]
        for v in deg2_vertices:
            if v not in g:
                continue
            u, w = g[v]
            if u != w:
                # Remove v
                del g[v]
                # Replace v with w in g[u]
                g[u] = [w if x == v else x for x in g[u]]
                # Replace v with u in g[w]
                g[w] = [u if x == v else x for x in g[w]]
                changed = True
                break

    # If the cleaned graph is empty or has no cycles, return empty set
    if not g or cyclomatic_number(g) == 0:
        return set()

    # Check for self-loops (which are cycles of length 1)
    for v in g:
        if v in g[v]:
            g_sub = {x: [y for y in neighbors if y != v] for x, neighbors in g.items() if x != v}
            sub_fvs = feedback_vertex_set_approx(g_sub, weights)
            sub_fvs.add(v)
            return sub_fvs

    # Check for multi-edges (which are cycles of length 2)
    for u in g:
        for v in g[u]:
            if g[u].count(v) > 1:
                eps = min(weights[u], weights[v])
                weights_next = dict(weights)
                weights_next[u] -= eps
                weights_next[v] -= eps
                
                g_sub = {x: list(neighbors) for x, neighbors in g.items()}
                sub_fvs = feedback_vertex_set_approx(g_sub, weights_next)
                
                # If neither is in sub_fvs, we must add the one that had weight reduced to 0
                if u not in sub_fvs and v not in sub_fvs:
                    if weights_next[u] <= 1e-9:
                        sub_fvs.add(u)
                    else:
                        sub_fvs.add(v)
                return sub_fvs

    # Normal local ratio step:
    # Compute cyclomatic weights delta_g(v) for all vertices in g
    cyc_g = cyclomatic_number(g)
    deltas = {}
    for v in g:
        g_without_v = {x: [y for y in neighbors if y != v] for x, neighbors in g.items() if x != v}
        deltas[v] = cyc_g - cyclomatic_number(g_without_v)
        if deltas[v] <= 0:
            deltas[v] = 1

    # Find epsilon = min_{v} weights[v] / deltas[v]
    eps = min(weights[v] / deltas[v] for v in g)
    
    # Subtract eps * deltas[v] from weights
    weights_next = dict(weights)
    for v in g:
        weights_next[v] -= eps * deltas[v]
        
    # Recurse
    sub_fvs = feedback_vertex_set_approx(g, weights_next)
    
    # Pruning / Extension Phase:
    zero_weight_vertices = [v for v in g if weights_next[v] <= 1e-9]
    for v in zero_weight_vertices:
        if v not in sub_fvs:
            sub_fvs.add(v)
            
    # Minimal FVS check
    final_fvs = set(sub_fvs)
    for v in list(final_fvs):
        if is_feedback_vertex_set(graph, final_fvs - {v}):
            final_fvs.remove(v)
            
    return final_fvs

def demo_feedback_vertex_set():
    print("=" * 60)
    print("Chapter 6: Feedback Vertex Set - 2-Approximation")
    print("=" * 60)
    
    # Example 1: Bipartite graph K_3,3 (6 vertices, 9 edges)
    # Vertices: 0,1,2 (left) and 3,4,5 (right)
    graph = {
        0: [3, 4, 5],
        1: [3, 4, 5],
        2: [3, 4, 5],
        3: [0, 1, 2],
        4: [0, 1, 2],
        5: [0, 1, 2]
    }
    weights = {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0}
    
    fvs = feedback_vertex_set_approx(graph, weights)
    print("\n1. Bipartite Graph K_3,3")
    print(f"  Graph vertices: {list(graph.keys())}")
    print(f"  Selected FVS: {fvs}")
    print(f"  FVS size: {len(fvs)} (Optimal is 2)")
    print(f"  Is valid FVS: {is_feedback_vertex_set(graph, fvs)}")

    # Example 2: Weighted graph
    # 0 - 1 - 2 - 0 (cycle 0-1-2-0)
    # 2 - 3 - 4 - 2 (cycle 2-3-4-2)
    # Vertex 2 is shared (weight 5.0), others have weight 2.0
    graph2 = {
        0: [1, 2],
        1: [0, 2],
        2: [0, 1, 3, 4],
        3: [2, 4],
        4: [2, 3]
    }
    weights2 = {0: 2.0, 1: 2.0, 2: 5.0, 3: 2.0, 4: 2.0}
    
    fvs2 = feedback_vertex_set_approx(graph2, weights2)
    print("\n2. Shared Center Cycle Graph")
    print(f"  Weights: {weights2}")
    print(f"  Selected FVS: {fvs2}")
    print(f"  Total weight: {sum(weights2[v] for v in fvs2)}")
    print(f"  Is valid FVS: {is_feedback_vertex_set(graph2, fvs2)}")

if __name__ == "__main__":
    demo_feedback_vertex_set()

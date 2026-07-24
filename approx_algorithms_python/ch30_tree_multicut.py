"""
Chapter 30: Multicut in Trees
=============================
Vazirani Ch. 18: Multicut and Integer Multicommodity Flow in Trees.
Implements:
1. Tree representation, path extraction, and LCA depth computation.
2. Primal-dual 2-approximation algorithm.
3. Reverse-delete pruning to guarantee the factor 2 bound.
"""

from typing import List, Tuple, Dict, Set

def get_tree_properties(
    n: int,
    edges: List[Tuple[int, int]],
    root: int = 0
) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Computes parent pointers and depths for all nodes in the tree."""
    adj = {i: [] for i in range(n)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        
    parent = {root: -1}
    depth = {root: 0}
    
    # BFS traversal to compute depths and parents
    queue = [root]
    while queue:
        curr = queue.pop(0)
        for neighbor in adj[curr]:
            if neighbor != parent[curr]:
                parent[neighbor] = curr
                depth[neighbor] = depth[curr] + 1
                queue.append(neighbor)
                
    return parent, depth

def get_lca(u: int, v: int, parent: Dict[int, int], depth: Dict[int, int]) -> int:
    """Finds the lowest common ancestor of u and v in the rooted tree."""
    while depth[u] > depth[v]:
        u = parent[u]
    while depth[v] > depth[u]:
        v = parent[v]
    while u != v:
        u = parent[u]
        v = parent[v]
    return u

def get_path_edges(u: int, v: int, parent: Dict[int, int]) -> List[Tuple[int, int]]:
    """Returns the list of edges on the unique path between u and v in the tree."""
    path_u = []
    curr = u
    while curr != -1:
        path_u.append(curr)
        curr = parent[curr]
        
    path_v = []
    curr = v
    while curr != -1:
        path_v.append(curr)
        curr = parent[curr]
        
    # Find LCA
    lca_node = -1
    set_v = set(path_v)
    for node in path_u:
        if node in set_v:
            lca_node = node
            break
            
    # Reconstruct path edges from u to lca, and lca to v
    edges_on_path = []
    curr = u
    while curr != lca_node:
        p = parent[curr]
        edges_on_path.append(tuple(sorted((curr, p))))
        curr = p
        
    curr = v
    while curr != lca_node:
        p = parent[curr]
        edges_on_path.append(tuple(sorted((curr, p))))
        curr = p
        
    return edges_on_path

def multicut_in_trees(
    n: int,
    edges: List[Tuple[int, int]],
    costs: List[float],
    pairs: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """
    Primal-dual 2-approximation for Multicut in Trees.
    """
    parent, depth = get_tree_properties(n, edges)
    
    # Sort pairs in descending order of the depth of their LCA
    pair_lcas = []
    for i, (u, v) in enumerate(pairs):
        lca = get_lca(u, v, parent, depth)
        pair_lcas.append((depth[lca], i, lca))
        
    pair_lcas.sort(key=lambda x: x[0], reverse=True)
    
    # Initialize dual variables and edge loads
    edge_to_idx = {tuple(sorted(e)): idx for idx, e in enumerate(edges)}
    load = [0.0] * len(edges)
    chosen_edges = []
    
    for _, pair_idx, lca in pair_lcas:
        u, v = pairs[pair_idx]
        path_edges = get_path_edges(u, v, parent)
        
        # Check if the path is already cut
        already_cut = False
        for e in path_edges:
            if e in chosen_edges:
                already_cut = True
                break
                
        if already_cut:
            continue
            
        # Grow dual y_i for this pair:
        # y_i = min_{e in path} (c_e - load(e))
        min_slack = float('inf')
        best_edge = None
        
        for e in path_edges:
            e_idx = edge_to_idx[e]
            slack = costs[e_idx] - load[e_idx]
            if slack < min_slack:
                min_slack = slack
                best_edge = e
            elif abs(slack - min_slack) < 1e-9:
                # Break tie by choosing the higher edge (closer to lca)
                # An edge (x, parent[x])'s height is determined by the depth of parent[x]
                u1, v1 = e
                u2, v2 = best_edge
                d1 = min(depth[u1], depth[v1])
                d2 = min(depth[u2], depth[v2])
                if d1 < d2:  # smaller depth means closer to root/lca
                    best_edge = e
                    
        # Update dual loads on path edges
        for e in path_edges:
            e_idx = edge_to_idx[e]
            load[e_idx] += min_slack
            
        if best_edge is not None:
            chosen_edges.append(best_edge)
            
    # Reverse-delete pruning
    pruned = list(chosen_edges)
    for e in reversed(chosen_edges):
        pruned.remove(e)
        
        # Check if all pairs are still cut by the pruned set
        all_cut = True
        for u, v in pairs:
            path_edges = get_path_edges(u, v, parent)
            cut_found = False
            for pe in path_edges:
                if pe in pruned:
                    cut_found = True
                    break
            if not cut_found:
                all_cut = False
                break
                
        if not all_cut:
            pruned.append(e)
            
    return pruned

def demo_tree_multicut():
    print("=" * 60)
    print("Chapter 30: Multicut in Trees")
    print("=" * 60)
    
    # Example tree
    # Root: 0
    # Edges: (0,1), (0,2), (1,3), (1,4), (2,5), (2,6)
    n = 7
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
    costs = [2.0, 4.0, 1.0, 3.0, 2.0, 2.0]
    
    # Pairs to disconnect
    pairs = [(3, 4), (5, 6), (3, 5)]
    
    print("\n1. Tree Instance and Demand Pairs:")
    print("  Edges & Costs:")
    for i, (u, v) in enumerate(edges):
        print(f"    Edge {i}: ({u}, {v}) cost={costs[i]}")
    print("  Pairs to disconnect:")
    for i, (u, v) in enumerate(pairs):
        print(f"    Pair {i}: ({u}, {v})")
        
    chosen = multicut_in_trees(n, edges, costs, pairs)
    print("\n2. Multicut Results:")
    print(f"  Selected Cut Edges: {chosen}")
    total_cost = sum(costs[edges.index(e) if e in edges else edges.index((e[1], e[0]))] for e in chosen)
    print(f"  Total Cut Cost:     {total_cost:.2f}")

if __name__ == "__main__":
    demo_tree_multicut()

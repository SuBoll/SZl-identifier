# MultiGraph implementation: supports parallel edges, no edges dropped
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['SimHei']      # or 'Microsoft YaHei'
matplotlib.rcParams['axes.unicode_minus'] = False        # display minus sign correctly

# -------------------------------
# Helper: print graph structure (sorted by idx)
# -------------------------------
def print_graph_structure(G):
    print("\n" + "="*50)
    print("Graph structure:")
    print("Vertices:", list(G.nodes()))

    print("Edges (with index):")
    edges = [(u, v, key, data) for u, v, key, data in G.edges(keys=True, data=True)]
    edges_sorted = sorted(edges, key=lambda t: t[3].get('idx', t[2]))
    for u, v, key, data in edges_sorted:
        print(f"  {u} -> {v}, key={key}, idx={data.get('idx', key)}")

    print(f"\nGraph info: n={G.number_of_nodes()}, m={G.number_of_edges()}, cycle space dim={G.number_of_edges()-G.number_of_nodes()+1}")
    print("="*50 + "\n")

# -------------------------------
# Visualization: MultiGraph with parallel edges
# -------------------------------
def visualize_multigraph(G: nx.MultiGraph, tree_edges, edge_vectors):
    fig, ax = plt.subplots(figsize=(8, 8), facecolor="#f8f9fa")
    ax.set_facecolor("#f8f9fa")

    # Layout: Petersen uses two concentric pentagons, others use kamada_kawai
    n, m = G.number_of_nodes(), G.number_of_edges()
    if n == 10 and m == 15:
        # Petersen graph: outer pentagon + inner pentagon
        angles = np.linspace(0, 2 * np.pi, 6)[:-1] - np.pi / 2
        pos = {i: (1.2 * np.cos(angles[i]), 1.2 * np.sin(angles[i])) for i in range(5)}
        pos.update({i + 5: (0.6 * np.cos(angles[i]), 0.6 * np.sin(angles[i])) for i in range(5)})
    else:
        pos = nx.kamada_kawai_layout(G)

    # Draw edges (arc offset for parallel edges, alternating sign)
    seen = {}
    for u, v, key, data in G.edges(keys=True, data=True):
        pair = (min(u, v), max(u, v))
        count = seen.get(pair, 0)
        seen[pair] = count + 1
        rad = 0.2 * (count - 1) * (-1) ** (count - 1) if count > 1 else 0
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=2.5,
                               connectionstyle=f"arc3, rad={rad}",
                               edge_color="#4a5568", alpha=0.85)

    # Draw nodes (with border for clarity)
    nx.draw_networkx_nodes(G, pos, node_size=600, node_color="#60a5fa",
                           edgecolors="#1e3a5f", linewidths=2)
    nx.draw_networkx_labels(G, pos, font_size=11, font_weight="bold",
                            font_color="#1e293b")

    ax.set_title("Graph structure", fontsize=14, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()
    plt.show()

# -------------------------------
# Format vector output (aligned, scientific notation supported)
# -------------------------------
def format_vector_aligned(vec, precision=4, width=12, sci=False):
    formatted = []
    for v in vec:
        if sci:
            formatted.append(f"{v:{width}.{precision}e}")
        else:
            formatted.append(f"{v:{width}.{precision}f}")
    return "[" + ", ".join(formatted) + "]"

# -------------------------------
# Core: assign cycle-space coefficient vectors to each edge (MultiGraph)
# -------------------------------
def assign_edge_vectors_multigraph(MG: nx.MultiGraph):
    n = MG.number_of_nodes()
    m = MG.number_of_edges()
    k = m - n + 1
    if k <= 0:
        raise ValueError("Cycle space dim k = m - n + 1 <= 0; graph has no cycles or is a tree.")

    all_edges = list(MG.edges(keys=True, data=True))
    edge_vectors = {(u, v, key): np.zeros(k, dtype=int) for u, v, key, _ in all_edges}

    # Index by endpoint pair: (min,max) -> [(u,v,key), ...] for tree edge matching
    pair_to_edges = {}
    for u, v, key, _ in all_edges:
        pair = (min(u, v), max(u, v))
        pair_to_edges.setdefault(pair, []).append((u, v, key))

    # Spanning tree and match MultiGraph edges
    T = nx.minimum_spanning_tree(MG)
    tree_edge_keys = []
    used_keys = set()
    for a, b in T.edges():
        pair = (min(a, b), max(a, b))
        candidates = pair_to_edges.get(pair, [])
        chosen = None
        for e in candidates:
            if e not in used_keys:
                chosen = e
                break
        if chosen is None and candidates:
            chosen = candidates[0]
        if chosen is None:
            raise RuntimeError(f"Cannot find MultiGraph edge for tree edge {(a,b)}.")
        tree_edge_keys.append(chosen)
        used_keys.add(chosen)

    # Non-tree edges
    non_tree_edges = [(u, v, key) for u, v, key, _ in all_edges if (u, v, key) not in used_keys]
    if len(non_tree_edges) < k:
        raise RuntimeError(f"Non-tree edge count {len(non_tree_edges)} < k={k}; check spanning tree logic.")

    non_tree_edges_sorted = sorted(non_tree_edges, key=lambda t: t[2])
    basis_edges = non_tree_edges_sorted[:k]
    for i, e in enumerate(basis_edges):
        edge_vectors[e][i] = 1

    # Incidence index: incident[v] = [(u,v,key), ...] for leaf peeling
    incident = {v: [] for v in MG.nodes()}
    for u, v, key in edge_vectors:
        incident[u].append((u, v, key))
        incident[v].append((u, v, key))

    # Tree edge vectors via leaf peeling
    unprocessed = set(tree_edge_keys)
    nodes = set(MG.nodes())
    while unprocessed:
        deg = {v: 0 for v in nodes}
        for u, v, key in unprocessed:
            deg[u] += 1
            deg[v] += 1
        leaves = [v for v in nodes if deg.get(v, 0) == 1]
        if not leaves:
            break
        for leaf in leaves:
            found_edge = None
            for u, v, key in list(unprocessed):
                if leaf in (u, v):
                    found_edge = (u, v, key)
                    break
            if found_edge is None:
                nodes.discard(leaf)
                continue
            u, v, key = found_edge
            child = leaf
            vec = np.zeros(k, dtype=int)
            for uu, vv, k2 in incident[child]:
                if uu == child:
                    vec += edge_vectors[(uu, vv, k2)]
                else:
                    vec -= edge_vectors[(uu, vv, k2)]
            edge_vectors[found_edge] = -vec if found_edge[0] == child else vec
            unprocessed.remove(found_edge)
            nodes.remove(child)

    return edge_vectors, tree_edge_keys

# -------------------------------
# Constrained optimization
# -------------------------------
def solve_with_constraints(MG: nx.MultiGraph, edge_vectors, d=2, p=2, repeats=20, maxiter=1000):
    edges_idx = sorted([(u, v, key, data) for u, v, key, data in MG.edges(keys=True, data=True)],
                       key=lambda t: t[3].get('idx', t[2]))
    edges_keys = [(u, v, key) for u, v, key, data in edges_idx]
    num_edges = len(edges_keys)
    k = len(next(iter(edge_vectors.values())))

    # Coefficient matrix C: (num_edges x k), each row = cycle-space coeffs for one edge; F = C @ X
    C = np.array([edge_vectors[e] for e in edges_keys], dtype=float)

    def make_cons(row_idx):
        row = C[row_idx]

        def cons(x):
            X = x.reshape((k, d))
            f = row @ X  # dot product
            return np.linalg.norm(f, ord=p) - 1.0
        return cons

    constraints = [{'type': 'ineq', 'fun': make_cons(i)} for i in range(num_edges)]

    best_M = np.inf
    best_x = None

    def objective(x):
        X = x.reshape((k, d))
        F = C @ X  # (num_edges x d)
        norms = np.linalg.norm(F, ord=p, axis=1)
        return float(np.max(norms))

    for t in range(repeats):
        x0=np.random.randn(k*d)
        try:
            res=minimize(objective,x0,method='SLSQP',constraints=constraints,
                         options={'maxiter':maxiter,'ftol':1e-6,'disp':False})
        except Exception as e:
            print(f"[Optimization error] Restart {t+1}: {e}")
            continue
        if res.success:
            M=float(res.fun)
            print(f"[Restart {t+1}/{repeats}] Success, M = {M:.6f}")
            if M<best_M:
                best_M=M
                best_x=res.x.copy()
        else:
            print(f"[Restart {t+1}/{repeats}] No convergence (message: {res.message})")

    if best_x is None:
        raise RuntimeError("No feasible solution found after all restarts")

    X = best_x.reshape((k, d))
    F = C @ X
    final_edge_vectors = {edges_keys[i]: F[i] for i in range(num_edges)}

    return best_M, final_edge_vectors

# -------------------------------
# Print edge vectors
# -------------------------------
def print_edge_vectors(MG, final_edge_vectors):
    print("\n" + "="*50)
    print("Edges (by idx):")
    print("-"*50)
    all_edges = sorted([(u, v, key, data) for u, v, key, data in MG.edges(keys=True, data=True)],
                      key=lambda t: t[3].get('idx', t[2]))
    for u, v, key, data in all_edges:
        vec = final_edge_vectors[(u, v, key)]
        norm = np.linalg.norm(vec)
        idx = data.get('idx', key)
        print(f" idx={idx:2d}  edge {u}--{v} (key={key}): vec = {format_vector_aligned(vec,precision=4,width=12,sci=True)}, ||·||_p = {norm:.6f}")
    print("="*50 + "\n")

# -------------------------------
# Print flow conservation
# -------------------------------
def print_flow_conservation(MG, final_edge_vectors, p=2):
    print("\n" + "="*50)
    print("Flow conservation (Out - In per node):")
    print("-"*50)
    d = len(next(iter(final_edge_vectors.values())))
    Out = {v: np.zeros(d) for v in MG.nodes()}
    In = {v: np.zeros(d) for v in MG.nodes()}
    for u, v, key, _ in MG.edges(keys=True, data=True):
        vec = final_edge_vectors[(u, v, key)]
        Out[u] += vec
        In[v] += vec
    for node in MG.nodes():
        diff = Out[node] - In[node]
        print(f" node {node}: diff = {format_vector_aligned(diff,precision=6,width=14,sci=True)}, ||diff||_p = {np.linalg.norm(diff,ord=p):.6e}")
    print("="*50 + "\n")

# -------------------------------
# Main example
# -------------------------------
def main_example():
    # Parameters: d=dim, p=norm, repeats=random restarts, maxiter=max iterations
    d = 2
    p = 2
    repeats = 100
    maxiter = 5000

    # Petersen graph
    n = 10
    edges = [(i, (i+1)%5) for i in range(5)] + [(i,i+5) for i in range(5)] + [(i+5,(i+2)%5+5) for i in range(5)]
    
    # Wheel graph W_k: n=k+1
    # n=4
    # edges=[(0,i) for i in range(1,n)] + [(i,i+1) for i in range(1,n-1)] + [(n-1,1)]

    # Cycle graph C_k
    # n = 5
    # edges = [(i,i+1) for i in range(n-1)] + [(4,0)]

    # Dipole graph Mu_k
    # n = 2
    # k = 5
    # edges = [(0,1) for i in range(k)]

    MG=nx.MultiGraph()
    MG.add_nodes_from(range(n))
    for idx,(u,v) in enumerate(edges):
        MG.add_edge(u,v,key=idx,idx=idx)

    print_graph_structure(MG)

    edge_vectors, tree_edge_keys=assign_edge_vectors_multigraph(MG)
    visualize_multigraph(MG, tree_edge_keys, edge_vectors)

    best_M, final_edge_vectors=solve_with_constraints(MG, edge_vectors,d=d,p=p,repeats=repeats,maxiter=maxiter)
    r=best_M+1.0

    print_edge_vectors(MG, final_edge_vectors)
    print_flow_conservation(MG, final_edge_vectors,p=p)
    print("\n" + "="*50)
    print("===== Final result =====")
    print(f"Min max edge norm M = {best_M:.6f}")
    print(f"Flow value r = {r:.6f}")
    print("="*50 + "\n")

if __name__=="__main__":
    main_example()

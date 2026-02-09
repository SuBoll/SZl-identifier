"""
Enumerate non-isomorphic n-vertex, m-edge multigraphs with constraints:
  - multiplicity per pair ≤ l-2
  - min degree ≥ l-1
  - edge connectivity ≥ 2
  - connected

Then test each for SZ_l using szl_solver (general l, odd or even).

Parameters n, m, l are configurable in main. Use m >= (n-1)(l-1).
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import sys
from typing import Dict, List, Tuple, Optional, TextIO

import importlib.util
import networkx as nx
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend; figures are saved only, not displayed
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np


# ---------- Dynamic import: SZlSolver (general modulus) ----------

def load_szl_solver_class() -> type:
    try:
        here = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        here = os.getcwd()
    solver_path = os.path.join(here, "szl_solver.py")
    spec = importlib.util.spec_from_file_location("szl_solver_module", solver_path)
    if spec is None or spec.loader is None:
        raise ImportError("Failed to load szl_solver.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # type: ignore[arg-type]
    spec.loader.exec_module(module)  # type: ignore[attr-defined]
    if not hasattr(module, "SZlSolver"):
        raise ImportError("SZlSolver not found in szl_solver.py")
    return module.SZlSolver  # type: ignore[return-value]


# ---------- Enumeration + isomorphism reduction ----------

def make_pairs(n: int) -> List[Tuple[int, int]]:
    """Return pairs (i,j) with 1<=i<j<=n in lexicographic order."""
    return [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]


def degree_from_counts(counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int) -> Dict[int, int]:
    deg: Dict[int, int] = {v: 0 for v in range(1, n + 1)}
    for (u, v), w in zip(pairs, counts):
        deg[u] += w
        deg[v] += w
    return deg


def is_connected_from_counts(counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int) -> bool:
    G = nx.Graph()
    G.add_nodes_from(range(1, n + 1))
    for (u, v), w in zip(pairs, counts):
        if w > 0:
            G.add_edge(u, v)
    return nx.is_connected(G)


def edge_connectivity_from_counts(counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int) -> int:
    Gm = nx.MultiGraph()
    Gm.add_nodes_from(range(1, n + 1))
    for (u, v), w in zip(pairs, counts):
        for _ in range(w):
            Gm.add_edge(u, v)
    return nx.edge_connectivity(Gm)


def pair_to_count(counts: Tuple[int, ...], pairs: List[Tuple[int, int]]) -> Dict[Tuple[int, int], int]:
    """Build dict (u,v) -> multiplicity for u < v."""
    return {p: c for p, c in zip(pairs, counts)}


def has_forbidden_3vertex_subgraph(
    counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int, l_val: int
) -> bool:
    """Check if graph contains a 3-vertex subgraph with ≥ 2l-2 edges (each pair mult ≤ l-1).
    Such 3-vertex graphs are known to be SZ_l; we exclude them to reduce enumeration.
    """
    if n < 3:
        return False
    w = pair_to_count(counts, pairs)
    threshold = 2 * l_val - 2
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            for k in range(j + 1, n + 1):
                m_ij = w.get((i, j), 0)
                m_ik = w.get((i, k), 0)
                m_jk = w.get((j, k), 0)
                if m_ij + m_ik + m_jk >= threshold:
                    return True
    return False


def _is_4vertex_szl_by_theorem(
    mult: Dict[Tuple[int, int], int], l_val: int
) -> bool:
    """Check if a 4-vertex graph (6 pairs) is SZ_l per the theorem.
    mult: (1,2),(1,3),(1,4),(2,3),(2,4),(3,4) -> multiplicity (already SZ_l-simplified).
    """
    # δ ≥ l-1
    deg = {
        1: mult.get((1, 2), 0) + mult.get((1, 3), 0) + mult.get((1, 4), 0),
        2: mult.get((1, 2), 0) + mult.get((2, 3), 0) + mult.get((2, 4), 0),
        3: mult.get((1, 3), 0) + mult.get((2, 3), 0) + mult.get((3, 4), 0),
        4: mult.get((1, 4), 0) + mult.get((2, 4), 0) + mult.get((3, 4), 0),
    }
    if min(deg.values()) < l_val - 1:
        return False
    # e ≥ 3l-3
    total = sum(mult.values())
    if total < 3 * l_val - 3:
        return False
    # (1) μ(G) = l-1
    mu = max(mult.values()) if mult else 0
    if mu == l_val - 1:
        return True
    # (2) ∃ 2-vertex set S with d(S) ≠ 2l-2
    four = [1, 2, 3, 4]
    for s in itertools.combinations(four, 2):
        sc = [x for x in four if x not in s]
        d_s = 0
        for u in s:
            for v in sc:
                a, b = (u, v) if u < v else (v, u)
                d_s += mult.get((a, b), 0)
        if d_s != 2 * l_val - 2:
            return True
    # (3) ∃ vertex v with d(v) ≢ l (mod 2)
    for v in four:
        if deg[v] % 2 != l_val % 2:
            return True
    return False


def has_forbidden_4vertex_subgraph(
    counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int, l_val: int
) -> bool:
    """Check if graph contains a 4-vertex induced subgraph that is SZ_l (by theorem).
    When n > 4, such graphs can be excluded (contain known SZ_l subgraph).
    """
    if n < 4:
        return False
    w = pair_to_count(counts, pairs)
    for quad in itertools.combinations(range(1, n + 1), 4):
        a, b, c, d = sorted(quad)
        mult_4 = {
            (1, 2): min(w.get((a, b), 0), l_val - 1),
            (1, 3): min(w.get((a, c), 0), l_val - 1),
            (1, 4): min(w.get((a, d), 0), l_val - 1),
            (2, 3): min(w.get((b, c), 0), l_val - 1),
            (2, 4): min(w.get((b, d), 0), l_val - 1),
            (3, 4): min(w.get((c, d), 0), l_val - 1),
        }
        if _is_4vertex_szl_by_theorem(mult_4, l_val):
            return True
    return False


def canonical_key(
    counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int
) -> Tuple[int, ...]:
    """Canonicalize an n-vertex multigraph by taking lexicographically smallest over all n! permutations."""
    w: Dict[Tuple[int, int], int] = {}
    for (u, v), val in zip(pairs, counts):
        a, b = (u, v) if u < v else (v, u)
        w[(a, b)] = val

    labels = list(range(1, n + 1))
    order_pairs = make_pairs(n)

    best: Optional[Tuple[int, ...]] = None
    for perm in itertools.permutations(labels):
        inv = {perm[i]: labels[i] for i in range(n)}
        tup = []
        for (a, b) in order_pairs:
            i_old = inv[a]
            j_old = inv[b]
            key = (i_old, j_old) if i_old < j_old else (j_old, i_old)
            tup.append(w[key])
        tup_t = tuple(tup)
        if best is None or tup_t < best:
            best = tup_t
    assert best is not None
    return best


def _enumerate_multiplicities_rec(
    pairs: List[Tuple[int, int]],
    remaining: int,
    max_per: int,
    prefix: List[int],
) -> List[Tuple[int, ...]]:
    """Recursively enumerate all multiplicity tuples summing to remaining, each ≤ max_per."""
    if len(pairs) == 0:
        if remaining == 0:
            return [tuple(prefix)]
        return []
    out: List[Tuple[int, ...]] = []
    hi = min(remaining, max_per)
    for k in range(hi, -1, -1):
        out.extend(
            _enumerate_multiplicities_rec(pairs[1:], remaining - k, max_per, prefix + [k])
        )
    return out


def enumerate_nonisomorphic_graphs(
    n: int, m: int, l_val: int
) -> List[Tuple[int, ...]]:
    """Enumerate n-vertex, m-edge non-isomorphic multigraphs.
    Constraints: mult ≤ l-2, min deg ≥ l-1, connected, edge connectivity ≥ 2,
    no 3-vertex subgraph with ≥ 2l-2 edges, and (when n>4) no 4-vertex SZ_l subgraph.
    """
    if l_val <= 0:
        raise ValueError("l must be positive.")
    min_degree = l_val - 1
    max_multiplicity = l_val - 2
    if max_multiplicity < 0:
        raise ValueError("l must be at least 2 (need l-2 >= 0 for max multiplicity).")

    pairs = make_pairs(n)

    all_tuples = _enumerate_multiplicities_rec(pairs, m, max_multiplicity, [])

    reps: Dict[Tuple[int, ...], Tuple[int, ...]] = {}
    for counts in all_tuples:
        deg = degree_from_counts(counts, pairs, n)
        if min(deg.values()) < min_degree:
            continue
        if not is_connected_from_counts(counts, pairs, n):
            continue
        if edge_connectivity_from_counts(counts, pairs, n) < 2:
            continue
        # Exclude graphs containing a 3-vertex subgraph with ≥ 2l-2 edges (known to be SZ_l).
        if has_forbidden_3vertex_subgraph(counts, pairs, n, l_val):
            continue
        # When n > 4, exclude graphs containing a 4-vertex SZ_l subgraph (by theorem).
        if n > 4 and has_forbidden_4vertex_subgraph(counts, pairs, n, l_val):
            continue
        key = canonical_key(counts, pairs, n)
        reps.setdefault(key, counts)

    return list(reps.values())


# ---------- Drawing ----------

def symmetric_rads(m: int, step: float) -> List[float]:
    if m <= 1:
        return [0.0]
    start = -step * (m - 1) / 2.0
    return [start + i * step for i in range(m)]


def draw_case(
    ax,
    counts: Tuple[int, ...],
    pairs: List[Tuple[int, int]],
    n: int,
    *,
    hub: int = 1,
    top: Optional[int] = None,
    radius: float = 0.8,
    node_size: int = 280,
    curve_step: float = 0.18,
    node_color: str = "#66D1B3",
    edge_color: str = "#5AA9E6",
) -> None:
    H = nx.Graph()
    H.add_nodes_from(range(1, n + 1))
    for (u, v), w in zip(pairs, counts):
        if w > 0:
            H.add_edge(u, v)

    pos: Dict[int, Tuple[float, float]] = {hub: (0.0, 0.0)}
    others = [v for v in range(1, n + 1) if v != hub]
    if top is not None and top in others:
        idx = others.index(top)
        others = others[idx:] + others[:idx]
    for i, v in enumerate(others):
        theta = math.pi / 2 + 2 * math.pi * i / max(1, len(others))
        pos[v] = (radius * math.cos(theta), radius * math.sin(theta))

    nx.draw_networkx_nodes(
        H, pos, node_color=node_color, node_size=node_size,
        linewidths=1.0, edgecolors="#2C3E50", ax=ax
    )
    nx.draw_networkx_labels(H, pos, font_color="#000000", font_size=9, ax=ax)

    for (u, v), mult in zip(pairs, counts):
        if mult == 0:
            continue
        rads = symmetric_rads(mult, curve_step)
        for rad in rads:
            r = rad
            if mult % 2 == 0 and abs(r) < 1e-6:
                r = 0.12
            patch = FancyArrowPatch(
                pos[u], pos[v],
                connectionstyle=f"arc3,rad={r}",
                arrowstyle="-",
                color=edge_color,
                linewidth=1.4,
                alpha=0.95,
                zorder=1,
            )
            ax.add_patch(patch)

    pad = 0.22
    lim = radius + pad
    ax.set_autoscale_on(False)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


# ---------- Build MultiGraph ----------

def build_multigraph(counts: Tuple[int, ...], pairs: List[Tuple[int, int]], n: int) -> nx.MultiGraph:
    Gm = nx.MultiGraph()
    Gm.add_nodes_from(range(1, n + 1))
    for (u, v), w in zip(pairs, counts):
        for _ in range(w):
            Gm.add_edge(u, v)
    return Gm


def log_print(msg: str, log_file: Optional[TextIO]) -> None:
    print(msg)
    if log_file is not None:
        log_file.write(msg + "\n")
        log_file.flush()


# ---------- Main ----------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Enumerate n-vertex, m-edge non-isomorphic multigraphs and test each for SZ_l (general l)."
    )
    parser.add_argument("--vertices", "-n", type=int, default=5, help="Number of vertices (default: 4).")
    parser.add_argument("--edges", "-m", type=int, default=21, help="Total edges (default: 12).")
    parser.add_argument("--modulus", "-l", type=int, default=6, help="Modulus l (default: 5).")
    parser.add_argument("--no-log", action="store_true", help="Do not write results to a txt file.")
    args = parser.parse_args()

    n = args.vertices
    m = args.edges
    l_val = args.modulus

    min_edges = (n - 1) * (l_val - 1)
    if m < min_edges:
        print("Edge count m={} is too small. Need m >= (n-1)(l-1) = {}.".format(m, min_edges))
        print("With such few edges, the graph cannot be SZ_l.")
        sys.exit(1)

    # Output folder: all output files go here
    output_dirname = "output_{}v{}e_l{}".format(n, m, l_val)
    output_dir = os.path.join(os.getcwd(), output_dirname)
    os.makedirs(output_dir, exist_ok=True)

    log_file: Optional[TextIO] = None
    if not args.no_log:
        results_filename = "nonisomorphic_{}v{}e_l{}_results.txt".format(n, m, l_val)
        results_path = os.path.join(output_dir, results_filename)
        try:
            log_file = open(results_path, "w", encoding="utf-8")
        except OSError:
            log_file = None

    def out(msg: str) -> None:
        log_print(msg, log_file)

    pairs = make_pairs(n)
    reps = enumerate_nonisomorphic_graphs(n, m, l_val)

    out("Output folder: {}".format(output_dir))
    out("Parameters: n={}, m={}, l={}".format(n, m, l_val))
    out("Constraints: min degree >= {}, max multiplicity per pair = {}, edge connectivity >= 2,".format(
        l_val - 1, l_val - 2))
    out("  no 3-vertex subgraph with >= 2l-2 edges (known SZ_l),")
    if n > 4:
        out("  no 4-vertex SZ_l subgraph (when n>4).")
    else:
        out("  (n=4: no 4-vertex filter; use full SZ_l test).")
    out("Number of non-isomorphic graphs: {}".format(len(reps)))

    if len(reps) == 0:
        out("No graphs satisfying the constraints.")
        if log_file is not None:
            log_file.close()
        return

    # ---------- Overview figures: 5x5 per image ----------
    graphs_per_fig = 25  # 5 * 5
    cols = 5
    n_figs = (len(reps) + graphs_per_fig - 1) // graphs_per_fig

    for fig_idx in range(n_figs):
        start = fig_idx * graphs_per_fig
        end = min(start + graphs_per_fig, len(reps))
        chunk = reps[start:end]
        n_chunk = len(chunk)

        fig_rows = max(1, (n_chunk + cols - 1) // cols)
        fig, axes = plt.subplots(fig_rows, cols, figsize=(cols * 2.2, fig_rows * 2.2))

        if isinstance(axes, np.ndarray):
            axes_flat = list(axes.flat)
        elif isinstance(axes, (list, tuple)):
            axes_flat = []
            for row in axes:
                if isinstance(row, (list, tuple, np.ndarray)):
                    axes_flat.extend(list(np.array(row).flat))
                else:
                    axes_flat.append(row)
        else:
            axes_flat = [axes]

        for i, counts in enumerate(chunk):
            ax = axes_flat[i]
            draw_case(ax, counts, pairs, n, hub=1, top=2 if n >= 2 else None, radius=0.8)
            ax.set_title("#{}".format(start + i + 1), fontsize=9)

        for j in range(n_chunk, len(axes_flat)):
            axes_flat[j].axis("off")

        plt.tight_layout()
        suffix = "_{}".format(fig_idx + 1) if n_figs > 1 else ""
        overview_filename = "nonisomorphic_{}v{}e_l{}_overview{}.png".format(n, m, l_val, suffix)
        overview_path = os.path.join(output_dir, overview_filename)
        fig.savefig(overview_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("Saved: {}".format(overview_path))  # console only, not written to txt

    # ---------- SZ_l testing ----------
    SZlSolver = load_szl_solver_class()
    non_sz_cases: List[Tuple[int, Tuple[int, ...], List[Dict[int, int]]]] = []

    for i, counts in enumerate(reps):
        Gm = build_multigraph(counts, pairs, n)
        solver = SZlSolver(Gm, l_val)
        infeasible_betas = solver.get_all_infeasible_betas(verbose=False)
        if infeasible_betas:
            non_sz_cases.append((i, counts, infeasible_betas))

    # Summary first: list which graphs are NOT SZ_l
    if non_sz_cases:
        indices = [c[0] + 1 for c in non_sz_cases]  # 1-based graph numbers
        out("")
        out("Graphs that are NOT SZ_{}: {}".format(l_val, ", ".join("#" + str(x) for x in indices)))
        out("")
        for i, counts, infeasible_betas in non_sz_cases:
            vertices = sorted(infeasible_betas[0].keys()) if infeasible_betas else []
            out("Graph #{} is NOT SZ_{}:".format(i + 1, l_val))
            out("- Edge multiplicities (sum={}): {}".format(m, dict(zip(pairs, counts))))
            deg = degree_from_counts(counts, pairs, n)
            out("- Vertex degrees: {}".format(deg))
            out("- All infeasible betas (vectors in vertex order {}):".format(vertices))
            for idx, beta in enumerate(infeasible_betas):
                vec = [beta[v] for v in vertices]
                out("  [{}] {}".format(idx + 1, vec))
            out("  Total: {} infeasible beta(s).".format(len(infeasible_betas)))
            out("")
    else:
        out("")
        out("All graphs are SZ_{}.".format(l_val))

    if log_file is not None:
        out("")
        out("Results written to: {}".format(results_path))
        log_file.close()


if __name__ == "__main__":
    main()

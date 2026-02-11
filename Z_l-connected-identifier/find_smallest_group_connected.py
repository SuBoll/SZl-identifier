"""
Find smallest simple 2-edge-connected graphs that are:
  (1) Z_2×Z_2-connected but NOT Z_4-connected,
  (2) Z_4-connected but NOT Z_2×Z_2-connected.

Only considers simple graphs (no multi-edges) and 2-edge-connected graphs.

Optimizations:
- Edge count: only consider graphs with m >= min{ceil(4n/3), ceil(3(n-1)/2)}.
- Superset (same vertex set): if G' = G - e (same n vertices, one edge removed)
  is 2-ec and A-connected, then G is A-connected; skip the solver.
"""

from __future__ import annotations

import itertools
import math
from typing import List, Set, Tuple

import networkx as nx

from group_connectivity_solver import GroupConnectivitySolver, _build_multigraph_from_edges


def min_edges_for_n(n: int) -> int:
    """Minimum edge count to consider: min{ceil(4n/3), ceil(3(n-1)/2)}."""
    return min(math.ceil(4 * n / 3), math.ceil(3 * (n - 1) / 2))


def all_pairs(n: int) -> List[Tuple[int, int]]:
    """List of all unordered pairs (i, j) with 1 <= i < j <= n, in fixed order."""
    return [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]


def enumerate_simple_2ec_graphs(n: int, m_min: int) -> List[List[Tuple[int, int]]]:
    """Enumerate simple 2-edge-connected graphs on vertices {1,...,n} with m >= m_min.

    Returns a list of edge lists. Each edge list is a list of (u, v) with u < v.
    """
    pairs = all_pairs(n)
    M = len(pairs)
    result: List[List[Tuple[int, int]]] = []
    # 2-ec requires at least n edges (e.g. cycle); we also require m >= m_min.
    for m in range(max(n, m_min), M + 1):
        for edge_indices in itertools.combinations(range(M), m):
            edges = [pairs[i] for i in edge_indices]
            Gm = _build_multigraph_from_edges(n, edges)
            G = nx.Graph(Gm)
            if not nx.is_connected(G):
                continue
            if nx.edge_connectivity(G) < 2:
                continue
            result.append(edges)
    return result


def edges_key(edges: List[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    """Canonical key for edge set (sorted tuple)."""
    return tuple(sorted(edges))


def is_2ec(n: int, edges: List[Tuple[int, int]]) -> bool:
    """True if the graph with given edges on n vertices is 2-edge-connected."""
    if len(edges) < n:
        return False
    Gm = _build_multigraph_from_edges(n, edges)
    return nx.is_connected(nx.Graph(Gm)) and nx.edge_connectivity(nx.Graph(Gm)) >= 2


def is_Z22_connected(edges: List[Tuple[int, int]], n: int) -> bool:
    """True if the simple graph with given edges (vertices 1..n) is Z_2×Z_2-connected."""
    Gm = _build_multigraph_from_edges(n, edges)
    solver = GroupConnectivitySolver(Gm, group_moduli=[2, 2])
    ok, _ = solver.is_A_connected(verbose=False)
    return ok


def is_Z4_connected(edges: List[Tuple[int, int]], n: int) -> bool:
    """True if the simple graph with given edges (vertices 1..n) is Z_4-connected."""
    Gm = _build_multigraph_from_edges(n, edges)
    solver = GroupConnectivitySolver(Gm, group_moduli=[4])
    ok, _ = solver.is_A_connected(verbose=False)
    return ok


def main() -> None:
    # 顶点数范围：只测试 n_min <= n <= n_max 的图（改这里即可）
    n_min = 6
    n_max = 10

    print("Enumerating simple 2-edge-connected graphs and testing Z_2×Z_2 vs Z_4 connectivity.")
    print(f"Vertex range: n from {n_min} to {n_max}.")
    print("Edge count: only m >= min{ceil(4n/3), ceil(3(n-1)/2)}.")
    print("Pruning: if a subgraph (remove one edge) is 2-ec and A-connected, skip solver.")
    print("Goal 1: smallest graph that is Z_2×Z_2-connected but NOT Z_4-connected.")
    print("Goal 2: smallest graph that is Z_4-connected but NOT Z_2×Z_2-connected.")
    print()

    found_z22_not_z4: List[Tuple[int, int, List[Tuple[int, int]]]] = []
    found_z4_not_z22: List[Tuple[int, int, List[Tuple[int, int]]]] = []
    stats_solver_z22 = 0
    stats_solver_z4 = 0
    stats_skip_z22 = 0
    stats_skip_z4 = 0

    for n in range(n_min, n_max + 1):
        m_min_n = min_edges_for_n(n)
        graphs = enumerate_simple_2ec_graphs(n, m_min_n)
        if not graphs:
            print(f"n={n}: no graphs with m>={m_min_n} (max edges = {n*(n-1)//2}).")
            continue

        by_edges: dict = {}
        for edges in graphs:
            m = len(edges)
            if m not in by_edges:
                by_edges[m] = []
            by_edges[m].append(edges)

        z22_set: Set[Tuple[Tuple[int, int], ...]] = set()
        z4_set: Set[Tuple[Tuple[int, int], ...]] = set()

        total = len(graphs)
        for idx, (m, edge_list) in enumerate(
            (m, edges) for m in sorted(by_edges.keys()) for edges in by_edges[m]
        ):
            edges = edge_list
            key = edges_key(edges)

            # Z_2×Z_2: 顶点集不变，若存在一条边 e 使得 G-e（同 n 点少一边）仍 2-ec 且已判为 Z_2×Z_2，则 G 必 Z_2×Z_2
            z22 = False
            for i in range(len(edges)):
                E_minus = edges[:i] + edges[i + 1 :]
                if is_2ec(n, E_minus) and edges_key(E_minus) in z22_set:
                    z22 = True
                    stats_skip_z22 += 1
                    break
            if not z22:
                z22 = is_Z22_connected(edges, n)
                stats_solver_z22 += 1
            if z22:
                z22_set.add(key)

            # Z_4: 同上（同顶点集，删一边）
            z4 = False
            for i in range(len(edges)):
                E_minus = edges[:i] + edges[i + 1 :]
                if is_2ec(n, E_minus) and edges_key(E_minus) in z4_set:
                    z4 = True
                    stats_skip_z4 += 1
                    break
            if not z4:
                z4 = is_Z4_connected(edges, n)
                stats_solver_z4 += 1
            if z4:
                z4_set.add(key)

            print(f"  n={n} m={m} graph {idx+1}/{total} Z22={z22} Z4={z4}", flush=True)
            if z22 and not z4:
                found_z22_not_z4.append((n, m, edges))
            if z4 and not z22:
                found_z4_not_z22.append((n, m, edges))

        print(f"n={n}: {len(graphs)} graphs (m>={m_min_n}); solver calls Z22={stats_solver_z22} Z4={stats_solver_z4}, skipped Z22={stats_skip_z22} Z4={stats_skip_z4}")
        if found_z22_not_z4:
            best = min(found_z22_not_z4, key=lambda x: (x[0], x[1]))
            print(f"  -> Found Z_2×Z_2 but not Z_4: n={best[0]}, m={best[1]}, edges={best[2]}")
        if found_z4_not_z22:
            best = min(found_z4_not_z22, key=lambda x: (x[0], x[1]))
            print(f"  -> Found Z_4 but not Z_2×Z_2: n={best[0]}, m={best[1]}, edges={best[2]}")

    print()
    print("========== Result ==========")
    if found_z22_not_z4:
        best = min(found_z22_not_z4, key=lambda x: (x[0], x[1]))
        print(f"Smallest Z_2×Z_2-connected but NOT Z_4-connected: n={best[0]}, m={best[1]}")
        print(f"  Edges: {best[2]}")
    else:
        print("No graph found that is Z_2×Z_2-connected but not Z_4-connected (in the range checked).")

    if found_z4_not_z22:
        best = min(found_z4_not_z22, key=lambda x: (x[0], x[1]))
        print(f"Smallest Z_4-connected but NOT Z_2×Z_2-connected: n={best[0]}, m={best[1]}")
        print(f"  Edges: {best[2]}")
    else:
        print("No graph found that is Z_4-connected but not Z_2×Z_2-connected (in the range checked).")


if __name__ == "__main__":
    main()

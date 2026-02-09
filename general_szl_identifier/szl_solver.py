"""
SZ_l decision and solving for general modulus l (≥1).

For general l: beta ∈ Z_{2l} with:
  - parity: beta(v) ≡ deg(v) (mod 2) for each vertex v
  - sum: sum_v beta(v) ≡ 0 (mod 2l)

A beta-orientation: for each vertex v, (outdegree - indegree) ≡ beta(v) (mod 2l).

SZ_l: for every legal beta, there exists a beta-orientation.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable

import networkx as nx


@dataclass(frozen=True)
class EdgeBundle:
    """Undirected vertex pair (no self-loops) with multiplicity, stored as (u, v, k) with u < v."""
    u: int
    v: int
    k: int


@dataclass
class OrientationSolution:
    """A beta-orientation solution for general l.

    - y_by_pair[(u,v)]: among k parallel edges between u and v (u < v), y edges oriented u -> v.
    - out_minus_in[v]: integer (outdegree - indegree) at vertex v.
    - beta: the given beta (values 0..2l-1).
    - directions: per-edge orientation list, each item is (tail, head).
    """
    modulus: int
    vertices: List[int]
    edge_bundles: List[EdgeBundle]
    y_by_pair: Dict[Tuple[int, int], int]
    out_minus_in: Dict[int, int]
    beta: Dict[int, int]
    directions: List[Tuple[int, int]]

    def pretty_print(self) -> None:
        mod2l = 2 * self.modulus
        print("—— One beta-orientation solution ——")
        print(f"l={self.modulus}, working mod 2l={mod2l}")
        print("beta (vector):", [self.beta[v] for v in self.vertices])
        print("Vertex out-in and check (mod 2l):")
        for v in self.vertices:
            val = self.out_minus_in[v]
            print(f"  v={v}: out-in={val}, mod {mod2l} = {val % mod2l}")
        print("Pairs (u,v), k, y(u->v), contribution to u (2y-k):")
        for eb in self.edge_bundles:
            y = self.y_by_pair[(eb.u, eb.v)]
            contrib = 2 * y - eb.k
            print(f"  ({eb.u},{eb.v}), k={eb.k}, y={y}, 2y-k={contrib}")
        print("Per-edge directions (first few):")
        for i, (a, b) in enumerate(self.directions[:40], 1):
            print(f"  e{i}: {a}->{b}")
        if len(self.directions) > 40:
            print(f"  ... total edges: {len(self.directions)}")


class SZlSolver:
    """SZ_l decision and solving for general modulus l (positive integer).

    Orientations: for pair {u,v} with multiplicity k, if y edges are oriented u -> v,
    then contribution at u is 2y - k, at v is -(2y - k). Contribution set per pair is
    {-k, -k+2, ..., k-2, k}.

    SZ_l: for every legal beta (parity + sum constraints over Z_{2l}),
    there exists a beta-orientation.
    """

    def __init__(self, multigraph: nx.MultiGraph, modulus: int):
        if modulus <= 0:
            raise ValueError("Modulus l must be a positive integer.")
        if any(u == v for u, v in multigraph.edges()):
            raise ValueError("Self-loops are not allowed.")
        if not nx.is_connected(nx.Graph(multigraph)):
            raise ValueError("The graph must be connected.")

        self.Gm: nx.MultiGraph = multigraph
        self.l: int = modulus
        self.mod2l: int = 2 * modulus

        self.vertices: List[int] = sorted(self.Gm.nodes())
        self.index_of_vertex: Dict[int, int] = {v: i for i, v in enumerate(self.vertices)}
        self.edge_bundles: List[EdgeBundle] = self._collect_edge_bundles()
        self.sign_by_vertex: List[List[Tuple[int, int]]] = self._build_signs()
        self.deg: Dict[int, int] = self._compute_degrees()

        # C_v = sum_e sign(v,e) * (-k_e), so that:
        #   C_v + 2 * sum_e sign(v,e) * y_e ≡ beta(v) (mod 2l)
        self.C_vec: List[int] = [
            sum(sign * (-self.edge_bundles[eidx].k) for eidx, sign in self.sign_by_vertex[v_idx])
            for v_idx in range(len(self.vertices))
        ]

    def _collect_edge_bundles(self) -> List[EdgeBundle]:
        count: Dict[Tuple[int, int], int] = {}
        for u, v in self.Gm.edges():
            a, b = (u, v) if u < v else (v, u)
            count[(a, b)] = count.get((a, b), 0) + 1
        return [EdgeBundle(u=a, v=b, k=k) for (a, b), k in sorted(count.items())]

    def _build_signs(self) -> List[List[Tuple[int, int]]]:
        n = len(self.vertices)
        sign_by_vertex: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
        for eidx, eb in enumerate(self.edge_bundles):
            u_idx = self.index_of_vertex[eb.u]
            v_idx = self.index_of_vertex[eb.v]
            sign_by_vertex[u_idx].append((eidx, +1))
            sign_by_vertex[v_idx].append((eidx, -1))
        return sign_by_vertex

    def _compute_degrees(self) -> Dict[int, int]:
        deg: Dict[int, int] = {v: 0 for v in self.vertices}
        for eb in self.edge_bundles:
            deg[eb.u] += eb.k
            deg[eb.v] += eb.k
        return deg

    # ---------- beta enumeration ----------

    def enumerate_betas(self) -> Iterable[Dict[int, int]]:
        """Enumerate all legal beta: V -> Z_{2l} with beta(v) ≡ deg(v) (mod 2) and sum(beta) ≡ 0 (mod 2l)."""
        n = len(self.vertices)
        mod = self.mod2l

        # Per vertex: choices are those x with x ≡ deg(v) (mod 2).
        choices_per_v: List[List[int]] = []
        for v in self.vertices:
            parity = self.deg[v] % 2
            choices = [x for x in range(mod) if x % 2 == parity]
            choices_per_v.append(choices)

        # Enumerate first n-1; last is determined by sum ≡ 0.
        for values in itertools.product(*[choices_per_v[i] for i in range(n - 1)]):
            s = sum(values) % mod
            last_parity = self.deg[self.vertices[-1]] % 2
            # last must satisfy: last ≡ s (mod 2) for sum parity, and last ≡ -s (mod 2l) for sum.
            # For sum ≡ 0: last ≡ -s (mod 2l). We need last ∈ [0..2l-1] with last ≡ last_parity (mod 2).
            last = (-s) % mod
            if last % 2 != last_parity:
                continue  # no valid last
            beta = {v: val for v, val in zip(self.vertices[:-1], values)}
            beta[self.vertices[-1]] = last
            yield beta

    # ---------- solving ----------

    def solve_for_beta(self, beta: Dict[int, int]) -> Tuple[bool, Optional[OrientationSolution]]:
        """Given legal beta, find a beta-orientation."""
        if set(beta.keys()) != set(self.vertices):
            raise ValueError("beta's vertex set does not match the graph's vertex set.")

        mod2l = self.mod2l
        mod_l = self.l

        # Legality checks
        for v in self.vertices:
            if beta[v] < 0 or beta[v] >= mod2l:
                return False, None
            if beta[v] % 2 != self.deg[v] % 2:
                return False, None
        if sum(beta.values()) % mod2l != 0:
            return False, None

        # gamma(v) = (beta(v) - C_v) / 2 in Z_l. Requires (beta - C_v) even (mod 2l).
        target_residue: List[int] = []
        for idx, v in enumerate(self.vertices):
            delta = (beta[v] - self.C_vec[idx]) % mod2l
            if delta % 2 != 0:
                return False, None
            target_residue.append((delta // 2) % mod_l)

        # Variables: y_e ∈ [0..k_e]. Constraint: sum_e sign(v,e)*y_e ≡ target_residue(v) (mod l).
        domains: List[range] = [range(0, eb.k + 1) for eb in self.edge_bundles]
        order = sorted(range(len(self.edge_bundles)), key=lambda eidx: len(domains[eidx]))

        y_sol: List[Optional[int]] = [None] * len(self.edge_bundles)
        n = len(self.vertices)
        partial_sum = [0] * n

        edge_u_idx: List[int] = [self.index_of_vertex[eb.u] for eb in self.edge_bundles]
        edge_v_idx: List[int] = [self.index_of_vertex[eb.v] for eb in self.edge_bundles]

        def remaining_range_for_vertex(vertex_idx: int, next_pos: int) -> Tuple[int, int]:
            L = 0
            U = 0
            for i_pos in range(next_pos, len(order)):
                eidx = order[i_pos]
                k = self.edge_bundles[eidx].k
                if edge_u_idx[eidx] == vertex_idx:
                    U += k
                elif edge_v_idx[eidx] == vertex_idx:
                    L -= k
            return L, U

        def ceil_div(a: int, b: int) -> int:
            return -((-a) // b)

        def residue_is_feasible(vertex_idx: int, next_pos: int) -> bool:
            L, U = remaining_range_for_vertex(vertex_idx, next_pos)
            need = (target_residue[vertex_idx] - (partial_sum[vertex_idx] % mod_l)) % mod_l
            qmin = ceil_div(L - need, mod_l)
            qmax = (U - need) // mod_l
            return qmin <= qmax

        def dfs(pos: int) -> bool:
            if pos == len(order):
                for v_idx in range(n):
                    if (partial_sum[v_idx] - target_residue[v_idx]) % mod_l != 0:
                        return False
                return True

            eidx = order[pos]
            eb = self.edge_bundles[eidx]
            u_idx = self.index_of_vertex[eb.u]
            v_idx = self.index_of_vertex[eb.v]

            for y in domains[eidx]:
                partial_sum[u_idx] += y
                partial_sum[v_idx] -= y

                pruned = False
                next_pos = pos + 1
                for vv in range(n):
                    if not residue_is_feasible(vv, next_pos):
                        pruned = True
                        break

                if not pruned and dfs(next_pos):
                    y_sol[eidx] = y
                    return True

                partial_sum[u_idx] -= y
                partial_sum[v_idx] += y

            return False

        ok = dfs(0)
        if not ok:
            return False, None

        # Assemble solution
        y_by_pair: Dict[Tuple[int, int], int] = {}
        out_minus_in: Dict[int, int] = {v: 0 for v in self.vertices}
        directions: List[Tuple[int, int]] = []

        for eidx, eb in enumerate(self.edge_bundles):
            y = y_sol[eidx]
            assert y is not None
            y_by_pair[(eb.u, eb.v)] = y
            contrib = 2 * y - eb.k
            out_minus_in[eb.u] += contrib
            out_minus_in[eb.v] -= contrib
            directions.extend([(eb.u, eb.v)] * y)
            directions.extend([(eb.v, eb.u)] * (eb.k - y))

        # Verify (mod 2l)
        for v in self.vertices:
            if out_minus_in[v] % mod2l != beta[v]:
                return False, None

        sol = OrientationSolution(
            modulus=self.l,
            vertices=self.vertices,
            edge_bundles=self.edge_bundles,
            y_by_pair=y_by_pair,
            out_minus_in=out_minus_in,
            beta=beta,
            directions=directions,
        )
        return True, sol

    def is_SZl(self, verbose: bool = False, max_beta: Optional[int] = None) -> Tuple[bool, Optional[Dict[int, int]]]:
        """Decide whether the graph is SZ_l. If not, return a witness infeasible beta."""
        cnt = 0
        for beta in self.enumerate_betas():
            cnt += 1
            ok, _ = self.solve_for_beta(beta)
            if verbose and cnt % 200 == 0:
                print(f"Checked betas: {cnt}")
            if not ok:
                return False, beta
            if max_beta is not None and cnt >= max_beta:
                break
        return True, None

    def get_all_infeasible_betas(
        self, verbose: bool = False, max_beta: Optional[int] = None
    ) -> List[Dict[int, int]]:
        """Enumerate all legal betas and return those for which no beta-orientation exists."""
        infeasible: List[Dict[int, int]] = []
        cnt = 0
        for beta in self.enumerate_betas():
            cnt += 1
            ok, _ = self.solve_for_beta(beta)
            if verbose and cnt % 200 == 0:
                print(f"Checked betas: {cnt}")
            if not ok:
                infeasible.append(beta)
            if max_beta is not None and cnt >= max_beta:
                break
        return infeasible


def build_graph_from_edges(n: int, edges: List[Tuple[int, int]]) -> nx.MultiGraph:
    Gm = nx.MultiGraph()
    Gm.add_nodes_from(list(range(1, n + 1)))
    for u, v in edges:
        Gm.add_edge(u, v)
    return Gm


def main():
    # ====== Example: edit the graph and l here (l can be any positive integer) ======
    mod = 5
    n = 4
    edges = (
        [(1, 2)] * 3
        + [(1, 3)] * 3
        + [(1, 4)] * 3
        + [(2, 3), (3, 4), (4, 2)]
    )
    # ================================================================================

    Gm = build_graph_from_edges(n, edges)
    solver = SZlSolver(Gm, mod)

    print("Graph info:")
    print(f"- vertices: {solver.vertices}")
    print(f"- total edges (with multiplicity): {sum(eb.k for eb in solver.edge_bundles)}")
    print(f"- degrees: {solver.deg}")
    print(f"- l={mod} (beta ∈ Z_{{2l}}, parity + sum constraints)")

    is_sz, witness = solver.is_SZl(verbose=True)
    if is_sz:
        print(f"\nConclusion: the graph IS SZ_{mod} ✓")
    else:
        print(f"\nConclusion: the graph is NOT SZ_{mod} ✗")
        print("One infeasible beta (vector):", [witness[v] for v in solver.vertices] if witness else None)

    # ====== Manual beta example ======
    # beta must be in 0..2l-1, beta(v) ≡ deg(v) mod 2, sum(beta) ≡ 0 mod 2l
    beta = {v: solver.deg[v] % 2 for v in solver.vertices}
    s = sum(beta.values()) % (2 * mod)
    if s != 0:
        beta = dict(beta)
        beta[solver.vertices[-1]] = (beta[solver.vertices[-1]] + 2 * mod - s) % (2 * mod)
    # ================================
    ok, sol = solver.solve_for_beta(beta)
    print("\nManual beta:", [beta[v] for v in solver.vertices])
    print("Feasibility:", "feasible" if ok else "infeasible")
    if ok and sol is not None:
        sol.pretty_print()


if __name__ == "__main__":
    main()

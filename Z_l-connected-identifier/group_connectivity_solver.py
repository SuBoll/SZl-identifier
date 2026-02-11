"""
Check whether a connected multigraph is A-connected for a finite abelian group A
(group connectivity: when a specific group A is given, we consider A-connectivity).

By the classification of finite abelian groups, every finite abelian group is
isomorphic to a direct product of cyclic groups. Hence we represent A as

    A ≅ Z_{m1} × Z_{m2} × ... × Z_{mr},

with each mi >= 2.  An element of A is then a tuple (x1, ..., xr) with
xi in {0, ..., mi-1}.  The group operation is component-wise addition mod mi.

Definitions (group-valued flows).
---------------------------------

Let A be a finite abelian group with identity 0_A, and let G be a connected
undirected multigraph (no self-loops).

- An A-boundary on G is a map

        beta : V(G) -> A

  such that the sum over all vertices is zero in A:

        sum_{v in V(G)} beta(v) = 0_A.

- Fix an arbitrary orientation of each undirected edge.  For an oriented edge
  e = (u -> v), a flow value is an element f(e) in A \ {0_A} (nowhere-zero).

  For a vertex v, define the net flow at v as

        net(v) = sum_{e outgoing from v} f(e)  -  sum_{e incoming to v} f(e)
               ∈ A,

  where subtraction is addition with the group inverse.

- A beta-flow is such an assignment f : E(G) -> A \ {0_A} satisfying

        net(v) = beta(v)   in A   for every vertex v.

  (Since A is abelian, the choice of the initial edge orientation does not
   matter: reversing an edge just replaces f(e) by -f(e).)

- The graph G is called A-connected if for every A-boundary beta there exists
  a beta-flow.

This module provides:

    GroupConnectivitySolver(Gm, group_moduli)

where

    - Gm is a networkx.MultiGraph (no self-loops, connected),
    - group_moduli is a list like [k]       for A = Z_k,
                           or [2, 2]   for A = Z_2 × Z_2,
                           or [m1,...,mr] in general.

The solver is exponential in both |V| and |E| (it enumerates all boundaries
and does a backtracking search for flows), and is meant for small graphs /
groups, similar in spirit to SZ_l experiments.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterable, Sequence

import networkx as nx


Elem = Tuple[int, ...]  # element of the finite abelian group


@dataclass(frozen=True)
class FiniteAbelianGroup:
    """Finite abelian group A ≅ Z_{m1} × ... × Z_{mr}, mi >= 2.

    Elements are tuples (x1,...,xr) with xi ∈ {0,...,mi-1}.
    The zero element is (0,...,0).
    """

    moduli: Tuple[int, ...]

    def __post_init__(self) -> None:
        if any(m <= 1 for m in self.moduli):
            raise ValueError("All moduli must be >= 2.")

    @property
    def rank(self) -> int:
        return len(self.moduli)

    def zero(self) -> Elem:
        return tuple(0 for _ in self.moduli)

    def add(self, a: Elem, b: Elem) -> Elem:
        return tuple((x + y) % m for (x, y, m) in zip(a, b, self.moduli))

    def neg(self, a: Elem) -> Elem:
        return tuple((-x) % m for (x, m) in zip(a, self.moduli))

    def elements(self) -> List[Elem]:
        """Return all elements of A as a list of tuples."""
        ranges: List[range] = [range(m) for m in self.moduli]
        return [tuple(coords) for coords in itertools.product(*ranges)]

    def is_zero(self, a: Elem) -> bool:
        return all(x == 0 for x in a)

    def format(self, a: Elem) -> str:
        """Pretty string for an element (scalar when rank=1)."""
        if self.rank == 1:
            return str(a[0])
        return "(" + ", ".join(str(x) for x in a) + ")"


@dataclass
class FlowSolution:
    """A beta-flow solution.

    - group: the underlying finite abelian group A.
    - vertices: list of vertex labels (sorted).
    - edges: list of oriented edges (tail, head); orientation is tail -> head.
    - beta: boundary map beta(v) in A.
    - flow: value f(e) in A for each edge index (nowhere-zero).
      flow[eidx] is the flow on the oriented edge edges[eidx] = (tail, head).
    """

    group: FiniteAbelianGroup
    vertices: List[int]
    edges: List[Tuple[int, int]]
    beta: Dict[int, Elem]
    flow: Dict[int, Elem]

    def get_orientation(self) -> List[Tuple[int, int]]:
        """Return the feasible orientation: list of (tail, head) for each edge."""
        return list(self.edges)

    def get_flow_values(self) -> Dict[int, Elem]:
        """Return flow value f(e) for each edge index (same keys as flow)."""
        return dict(self.flow)

    def pretty_print(self, max_edges: Optional[int] = None) -> None:
        """Print boundary, feasible orientation, and flow values.

        If max_edges is set, only the first max_edges edges are printed (default: all).
        """
        print("—— beta-flow 解：定向与流值 ——")
        print(f"A ≅ " + " × ".join(f"Z_{m}" for m in self.group.moduli))
        print("边界 beta (按顶点):", [self.group.format(self.beta[v]) for v in self.vertices])
        print()
        print("定向 (orientation)：每条边为 tail -> head")
        n_show = len(self.edges) if max_edges is None else min(max_edges, len(self.edges))
        for eidx in range(n_show):
            u, v = self.edges[eidx]
            print(f"  边 {eidx}: {u} -> {v}")
        if len(self.edges) > n_show:
            print(f"  ... 共 {len(self.edges)} 条边")
        print()
        print("流值 (flow)：f(e) ∈ A \\ {{0}}")
        for eidx in range(n_show):
            u, v = self.edges[eidx]
            val = self.flow[eidx]
            print(f"  边 {eidx} ({u}->{v}): f = {self.group.format(val)}")
        if len(self.edges) > n_show:
            print(f"  ... 共 {len(self.edges)} 条边")


class GroupConnectivitySolver:
    """Check A-connectivity of a multigraph for a finite abelian group A
    (group connectivity: given group A, we consider A-connectivity).

    Parameters
    ----------
    multigraph : nx.MultiGraph
        Connected multigraph without self-loops.
    group_moduli : Sequence[int]
        List of moduli [m1, ..., mr] describing the abelian group
        A ≅ Z_{m1} × ... × Z_{mr}, with each mi >= 2.

    Notes
    -----
    - Complexity is exponential in both |V| and |E|.
    - Intended for small graphs and small groups (e.g., experiments).
    """

    def __init__(self, multigraph: nx.MultiGraph, group_moduli: Sequence[int]):
        if any(u == v for u, v in multigraph.edges()):
            raise ValueError("Self-loops are not allowed.")
        if not nx.is_connected(nx.Graph(multigraph)):
            raise ValueError("The graph must be connected.")
        if not group_moduli:
            raise ValueError("group_moduli must be a non-empty list of integers >= 2.")

        self.Gm: nx.MultiGraph = multigraph
        self.group = FiniteAbelianGroup(tuple(int(m) for m in group_moduli))

        # Vertex ordering and index map.
        self.vertices: List[int] = sorted(self.Gm.nodes())
        self.index_of_vertex: Dict[int, int] = {v: i for i, v in enumerate(self.vertices)}
        self.n_vertices: int = len(self.vertices)

        # Fix an arbitrary orientation for each parallel edge: tail < head.
        # Each parallel edge is kept as a separate variable; we do not bundle.
        self.edges: List[Tuple[int, int]] = []
        for u, v in self.Gm.edges():
            a, b = (u, v) if u < v else (v, u)
            self.edges.append((a, b))
        self.m_edges: int = len(self.edges)

        # For each vertex v_idx, store list of (edge_index, sign) where
        # sign = +1 if v is the tail of that edge, -1 if v is the head.
        self.sign_by_vertex: List[List[Tuple[int, int]]] = [[] for _ in range(self.n_vertices)]
        for eidx, (u, v) in enumerate(self.edges):
            u_idx = self.index_of_vertex[u]
            v_idx = self.index_of_vertex[v]
            self.sign_by_vertex[u_idx].append((eidx, +1))
            self.sign_by_vertex[v_idx].append((eidx, -1))

        # Precompute all group elements and non-zero elements.
        self._all_elems: List[Elem] = self.group.elements()
        self._zero: Elem = self.group.zero()
        self._nonzero_elems: List[Elem] = [a for a in self._all_elems if not self.group.is_zero(a)]

        if not self._nonzero_elems:
            raise ValueError("Group has no non-zero elements (this should not happen).")

    # ------------------------------------------------------------------
    #  A-boundaries
    # ------------------------------------------------------------------

    def enumerate_boundaries(self) -> Iterable[Dict[int, Elem]]:
        """Enumerate all A-boundaries beta : V -> A with sum_v beta(v) = 0_A.

        The enumeration is done by freely choosing values for the first n-1
        vertices, and setting the last vertex to make the total sum zero.
        """
        n = self.n_vertices
        verts = self.vertices
        elems = self._all_elems
        zero = self._zero

        if n == 0:
            return  # empty graph not expected here

        # If there is only one vertex, we must have beta(v) = 0 for the sum to be zero.
        if n == 1:
            yield {verts[0]: zero}
            return

        # Choose beta for vertices v0..v_{n-2} freely, then fix v_{n-1}.
        first_vertices = verts[:-1]
        last_vertex = verts[-1]

        for values in itertools.product(elems, repeat=n - 1):
            total = zero
            for val in values:
                total = self.group.add(total, val)
            # beta(last) must satisfy total + beta(last) = 0  =>  beta(last) = -total.
            beta_last = self.group.neg(total)

            beta: Dict[int, Elem] = {
                v: val for v, val in zip(first_vertices, values)
            }
            beta[last_vertex] = beta_last
            yield beta

    # ------------------------------------------------------------------
    #  Flow solving for a fixed boundary
    # ------------------------------------------------------------------

    def _check_beta_is_boundary(self, beta: Dict[int, Elem]) -> None:
        """Validate that beta is a proper A-boundary on this graph.

        Conditions:
        - beta is defined exactly on self.vertices;
        - sum_v beta(v) = 0 in A.
        """
        if set(beta.keys()) != set(self.vertices):
            raise ValueError("beta must be defined on exactly the vertex set of the graph.")

        total = self._zero
        for v in self.vertices:
            total = self.group.add(total, beta[v])
        if not self.group.is_zero(total):
            raise ValueError("beta is not an A-boundary: sum_v beta(v) != 0 in A.")

    def solve_for_boundary(self, beta: Dict[int, Elem]) -> Optional[FlowSolution]:
        """Try to find a nowhere-zero beta-flow for the given boundary beta.

        Returns a FlowSolution if successful, or None if no flow exists.
        """
        self._check_beta_is_boundary(beta)

        n = self.n_vertices

        # partial_sum[v_idx] is the current net(v) = sum_out f - sum_in f in A,
        # considering only edges assigned so far.
        partial_sum: List[Elem] = [self._zero] * n
        flow_values: List[Optional[Elem]] = [None] * self.m_edges

        # Order edges by min(deg(u), deg(v)) ascending to encourage pruning early.
        edge_order = sorted(
            range(self.m_edges),
            key=lambda eidx: min(
                len(self.sign_by_vertex[self.index_of_vertex[self.edges[eidx][0]]]),
                len(self.sign_by_vertex[self.index_of_vertex[self.edges[eidx][1]]]),
            ),
        )
        # For pruning: if a vertex has all its incident edges in edge_order[0:pos] and partial_sum != beta, fail.
        incident_edge_indices: List[List[int]] = [
            [eidx for (eidx, _) in self.sign_by_vertex[v_idx]]
            for v_idx in range(n)
        ]

        def dfs(pos: int) -> bool:
            if pos == self.m_edges:
                for v_idx, v in enumerate(self.vertices):
                    if partial_sum[v_idx] != beta[v]:
                        return False
                return True

            eidx = edge_order[pos]
            u, v = self.edges[eidx]
            u_idx = self.index_of_vertex[u]
            v_idx = self.index_of_vertex[v]

            for val in self._nonzero_elems:
                old_u = partial_sum[u_idx]
                old_v = partial_sum[v_idx]

                partial_sum[u_idx] = self.group.add(partial_sum[u_idx], val)
                partial_sum[v_idx] = self.group.add(partial_sum[v_idx], self.group.neg(val))

                flow_values[eidx] = val

                # Prune: any vertex whose incident edges are all in edge_order[0:pos+1] must have partial_sum == beta.
                assigned = set(edge_order[: pos + 1])
                ok_so_far = True
                for vv_idx in range(n):
                    if not all(eid in assigned for eid in incident_edge_indices[vv_idx]):
                        continue
                    if partial_sum[vv_idx] != beta[self.vertices[vv_idx]]:
                        ok_so_far = False
                        break
                if ok_so_far and dfs(pos + 1):
                    return True

                partial_sum[u_idx] = old_u
                partial_sum[v_idx] = old_v
                flow_values[eidx] = None

            return False

        ok = dfs(0)
        if not ok:
            return None

        flow_dict: Dict[int, Elem] = {
            eidx: flow_values[eidx] for eidx in range(self.m_edges)  # type: ignore[arg-type]
        }

        return FlowSolution(
            group=self.group,
            vertices=self.vertices,
            edges=self.edges,
            beta=beta,
            flow=flow_dict,
        )

    # ------------------------------------------------------------------
    #  A-connectivity tests
    # ------------------------------------------------------------------

    def is_A_connected(self, verbose: bool = False) -> Tuple[bool, Optional[Dict[int, Elem]]]:
        """Check whether the graph is A-connected.

        Returns
        -------
        (True, None)
            if for every A-boundary beta there exists a nowhere-zero beta-flow.

        (False, beta)
            if a counterexample beta was found (no beta-flow exists for this beta).
        """
        for beta in self.enumerate_boundaries():
            sol = self.solve_for_boundary(beta)
            if sol is None:
                if verbose:
                    print("Found boundary with no beta-flow:")
                    print({v: self.group.format(beta[v]) for v in self.vertices})
                return False, beta
        return True, None

    def get_all_infeasible_boundaries(self, verbose: bool = False) -> List[Dict[int, Elem]]:
        """Return a list of all A-boundaries beta with no beta-flow.

        Warning: potentially enormous; use only for very small graphs/groups.
        """
        infeasible: List[Dict[int, Elem]] = []
        for beta in self.enumerate_boundaries():
            sol = self.solve_for_boundary(beta)
            if sol is None:
                infeasible.append(beta)
                if verbose:
                    print("Infeasible boundary beta:")
                    print({v: self.group.format(beta[v]) for v in self.vertices})
        return infeasible

    def check_beta(
        self,
        beta: Dict[int, Elem],
        verbose: bool = True,
        max_edges_in_print: Optional[int] = None,
    ) -> Tuple[bool, Optional[FlowSolution]]:
        """Check whether the given beta admits a nowhere-zero beta-flow.

        If a flow exists and verbose is True, prints the feasible orientation
        and flow values.

        Parameters
        ----------
        beta : dict
            A-boundary (sum of beta(v) must be 0_A).
        verbose : bool
            If True and feasible, print orientation and flow; if infeasible, print a message.
        max_edges_in_print : int, optional
            If set, only the first this many edges are printed (default: all).

        Returns
        -------
        (True, solution)
            If a beta-flow exists; solution contains .get_orientation() and .get_flow_values().
        (False, None)
            If no beta-flow exists for this beta.
        """
        sol = self.solve_for_boundary(beta)
        if sol is None:
            if verbose:
                print("该 beta 不可行：不存在 nowhere-zero beta-flow。")
                print("beta =", {v: self.group.format(beta[v]) for v in self.vertices})
            return False, None
        if verbose:
            sol.pretty_print(max_edges=max_edges_in_print)
        return True, sol


def _build_multigraph_from_edges(n: int, edges: Sequence[Tuple[int, int]]) -> nx.MultiGraph:
    """Helper to build a MultiGraph on vertices {1,...,n} from a list of edges."""
    Gm = nx.MultiGraph()
    for v in range(1, n + 1):
        Gm.add_node(v)
    for u, v in edges:
        if u == v:
            raise ValueError("Self-loops are not allowed.")
        Gm.add_edge(u, v)
    return Gm


def main_demo() -> None:
    """Small demo: check A-connectivity and, for a given beta, output orientation and flow."""
    # Example graph: a triangle 1-2-3-1.
    n = 3
    edges = [(1, 2), (2, 3), (3, 1)]
    Gm = _build_multigraph_from_edges(n, edges)

    # Example 1: A = Z_3.
    print("=== Demo: A = Z_3 ===")
    solver_z3 = GroupConnectivitySolver(Gm, group_moduli=[3])
    ok, bad_beta = solver_z3.is_A_connected(verbose=True)
    print("Is A-connected (A = Z_3)?", ok)

    # Check a specific beta and, if feasible, output orientation and flow.
    print("\n--- 对指定 beta 检查并输出定向与流值 ---")
    beta_example = {1: (1,), 2: (1,), 3: (1,)}  # 1+1+1 = 3 ≡ 0 (mod 3)
    feasible, sol = solver_z3.check_beta(beta_example, verbose=True)
    if feasible and sol is not None:
        orient = sol.get_orientation()
        flow_vals = sol.get_flow_values()
        print("定向 (list):", orient)
        print("流值 (dict eidx -> value):", {e: solver_z3.group.format(flow_vals[e]) for e in flow_vals})

    # Example 2: A = Z_2 × Z_2.
    print("\n=== Demo: A = Z_2 × Z_2 ===")
    solver_z2z2 = GroupConnectivitySolver(Gm, group_moduli=[2, 2])
    ok2, bad_beta2 = solver_z2z2.is_A_connected(verbose=True)
    print("Is A-connected (A = Z_2 × Z_2)?", ok2)
    if not ok2:
        print("Counterexample beta:", {v: solver_z2z2.group.format(bad_beta2[v]) for v in solver_z2z2.vertices})

    # 对 Z_2 × Z_2 给一个 beta 例子并输出对应的流
    print("\n--- Z_2×Z_2：对指定 beta 检查并输出定向与流值 ---")
    beta_z2z2 = {1: (1, 0), 2: (0, 1), 3: (1, 1)}  # (1,0)+(0,1)+(1,1) = (0,0) in Z_2×Z_2
    feasible2, sol2 = solver_z2z2.check_beta(beta_z2z2, verbose=True)
    if feasible2 and sol2 is not None:
        print("定向 (list):", sol2.get_orientation())
        print("流值 (dict eidx -> value):", {e: solver_z2z2.group.format(sol2.get_flow_values()[e]) for e in sol2.get_flow_values()})


if __name__ == "__main__":
    main_demo()

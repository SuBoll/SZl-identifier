## $SZ_l$ (general modulus) experiments: solver + enumerator

This folder contains scripts for **general modulus $l$** (odd or even), complementary to the odd-modulus version (`szl_odd_solver.py` + `generate_nonisomorphic_4v_szl.py`):

- `szl_solver.py`: an $SZ_l$ decision and $\beta$-orientation solver for **arbitrary positive** modulus $l$.
- `generate_nonisomorphic_szl.py`: for **given $n$, $m$, $l$**, enumerates $n$-vertex, $m$-edge (with multiplicity) undirected multigraphs under constraints (min degree $\ge l-1$, max multiplicity per pair $\le l-2$, edge connectivity $\ge 2$, connected); excludes graphs containing known $SZ_l$ subgraphs (3-vertex or 4-vertex); removes isomorphic duplicates; renders overview figures; and tests each graph for $SZ_l$; for graphs that are not $SZ_l$, it prints **all** infeasible $\beta$.

This README is the **main document**: it explains how the two scripts fit together, then gives a **step-by-step** account of the mathematical modeling and implementation workflow (aligned with the actual code), so results can be reviewed and reproduced.

> Note: you may see `libpng warning: iCCP: known incorrect sRGB profile` at runtime. This is typically a benign PNG color-profile warning and does not affect saving or displaying the figure.

---

## 0. Recommended layout

Keep these files in the same directory:

- `szl_solver.py`
- `generate_nonisomorphic_szl.py`
- `README_CN.md` (Chinese)
- `README.md` (this file, English)

`generate_nonisomorphic_szl.py` dynamically imports `szl_solver.py` from the same directory to access `SZlSolver`, so they must be co-located unless you adjust the import path.

Running the script creates an `output_{n}v{m}e_l{l}` folder in the current working directory, containing overview PNGs and results TXT.

---

## 1. Requirements

- Python 3.8+
- `networkx`
- `matplotlib`
- `numpy`

Install example:

```bash
python -m pip install networkx matplotlib numpy
```

---

## 2. Quickstart (recommended)

### 2.1 Run enumeration + overview figures + $SZ_l$ testing

```bash
python generate_nonisomorphic_szl.py [--vertices n] [--edges m] [--modulus l]
```

Defaults are `n=5`, `m=21`, `l=6`. Example: `python generate_nonisomorphic_szl.py -n 4 -m 12 -l 5`.

This will:

- create an `output_{n}v{m}e_l{l}` folder in the current working directory
- save overview PNG(s) inside the folder (up to 5×5=25 graphs per image; multiple images if needed)
- print $SZ_l$ testing results to the console; if some graphs are not $SZ_l$, it first prints a summary of their indices, then prints **all** infeasible $\beta$ (as vectors) for each
- write the same run output (excluding "Saved: ..." paths) to a results file `nonisomorphic_{n}v{m}e_l{l}_results.txt` in the same folder
- use `--no-log` to disable writing the results file

Note: figures are saved only; no window is displayed (Agg backend).

### 2.2 Run the solver demo only

```bash
python szl_solver.py
```

The script contains a built-in `main()` demo. You can edit `mod/n/edges/beta` to test different graphs and boundaries.

---

## 3. `szl_solver.py`: problem statement

### 3.1 Valid boundaries $\beta$ (general modulus $l$)

Let $G$ be a connected undirected multigraph (parallel edges allowed, no self-loops) and let $l$ be a positive integer.

A function $\beta:V(G)\to\mathbb{Z}_{2l}$ is a valid boundary if:

- $\beta(v)\in\{0,1,\dots,2l-1\}$ for all $v$
- $\beta(v)\equiv \deg(v)\pmod 2$ (parity constraint)
- $\sum_{v\in V}\beta(v)\equiv 0\pmod{2l}$

The parity constraint ensures that $\gamma(v)$ in the modeling is well-defined (division by 2); the sum condition is necessary since $\sum_v(\text{out}-\text{in})=0$ for any orientation.

### 3.2 $\beta$-orientation and $SZ_l$

A $\beta$-orientation is an orientation of all edges such that for each vertex $v$:

$$
(\mathrm{outdeg}(v)-\mathrm{indeg}(v))\equiv \beta(v)\pmod{2l}
$$

The graph $G$ is $SZ_l$ if such an orientation exists for **every** valid boundary $\beta$.

---

## 4. `szl_solver.py`: modeling (as implemented)

### 4.1 Bundling parallel edges: one integer $y$ per vertex pair

Same as the odd-modulus case: for an unordered pair $\{u,v\}$ (normalized as $u<v$) with multiplicity $k$, define $y\in\{0,1,\dots,k\}$ as the number of edges oriented $u\to v$. Then:

- net contribution at $u$ (out-in) is $2y-k$
- net contribution at $v$ is $-(2y-k)$

### 4.2 Signs and constant terms $C_v$

Same as the odd-modulus case:

$$
C_v+2S(v)\equiv \beta(v)\pmod{2l},\qquad
C_v=\sum_{e\ni v}\mathrm{sign}(v,e)\cdot(-k_e)
$$

### 4.3 General modulus: dividing by 2 to obtain $\gamma$

When $l$ is even, 2 is not invertible in $\mathbb{Z}_{2l}$. The parity constraint $\beta(v)\equiv \deg(v)\pmod 2$ ensures that $(\beta(v)-C_v)$ is even (mod $2l$). Define:

$$
\gamma(v)=\frac{(\beta(v)-C_v)\bmod 2l}{2}\in\mathbb{Z}_l
$$

The constraint becomes:

$$
S(v)\equiv \gamma(v)\pmod l
$$

In the code, the right-hand side is precomputed as `target_residue[v_idx]` (aligned with the sorted vertex list `vertices`).

---

## 5. `szl_solver.py`: implementation workflow (step-by-step)

### 5.1 Core data structures

- `EdgeBundle(u, v, k)`: a normalized vertex pair with multiplicity $k$ (always `u < v`).
- `OrientationSolution`: stores a solution, including `y_by_pair`, `out_minus_in`, `directions` (verification is mod $2l$).

### 5.2 `SZlSolver.__init__`: preprocessing (step-by-step)

1. **Input validation**: `modulus > 0`; no self-loops; graph connected (checked via `nx.is_connected(nx.Graph(multigraph))`).
2. **Vertices and indices**: `vertices = sorted(self.Gm.nodes())`, `index_of_vertex[v] = i` (sorted by vertex label).
3. **`_collect_edge_bundles()`**:
   - Iterate `self.Gm.edges()`; for each edge `(u,v)` normalize to `(a,b)` with `a,b = (u,v) if u<v else (v,u)`;
   - Accumulate multiplicity in `count[(a,b)]`;
   - Return `[EdgeBundle(u=a, v=b, k=k) for (a,b), k in sorted(count.items())]`, so pairs are in lexicographic order.
4. **`_build_signs()`**:
   - For each bundle `e=(u,v)` (with u<v), set `sign(u,e)=+1`, `sign(v,e)=-1`;
   - For each vertex `v_idx`, build `sign_by_vertex[v_idx] = [(eidx, +1 or -1), ...]` listing incident bundles and their signs.
5. **`_compute_degrees()`**: For each bundle's `k`, add to `deg[u]` and `deg[v]`.
6. **`C_vec`**: For vertex `v_idx`, `C_vec[v_idx] = sum(sign * (-k_e) for eidx, sign in sign_by_vertex[v_idx])`, i.e., the precomputed $C_v$.

### 5.3 `enumerate_betas()`: enumerate all valid $\beta$ (step-by-step)

1. For each vertex `v`, admissible values are $\{x : x\equiv \deg(v)\pmod 2,\; x\in [0,2l)\}$, i.e. `choices_per_v[i] = [x for x in range(mod2l) if x % 2 == deg[v] % 2]`.
2. Use `itertools.product(*[choices_per_v[i] for i in range(n-1)])` to enumerate the first $n-1$ vertices.
3. For each `values`, compute `s = sum(values) % mod2l` and set `last = (-s) % mod2l` so the total is 0 (mod $2l$).
4. Check `last % 2 == deg[vertices[-1]] % 2`; if not equal, skip (no valid last).
5. Build `beta = {v: val for v, val in zip(vertices[:-1], values)}`, set `beta[vertices[-1]] = last`, and `yield beta`.

### 5.4 `solve_for_beta(beta)`: DFS + pruning (step-by-step)

1. **Validity check**: vertex set matches; for each `v`, `0 <= beta[v] < 2l` and `beta[v] % 2 == deg[v] % 2`; `sum(beta.values()) % mod2l == 0`.
2. **Compute target_residue (i.e. $\gamma$)**: For each vertex `v`, `delta = (beta[v] - C_vec[idx]) % mod2l`; if `delta % 2 != 0`, return infeasible immediately; else `target_residue[idx] = (delta // 2) % l`.
3. **Variable domains**: `domains[eidx] = range(0, k_e + 1)`, i.e. $y_e \in \{0,\dots,k_e\}$.
4. **Variable ordering**: `order = sorted(range(len(edge_bundles)), key=lambda eidx: len(domains[eidx]))`, so smaller-domain bundles are assigned first (better pruning).
5. **DFS state**: `partial_sum[v_idx]` is the integer sum of $S(v)$ from assigned bundles; `y_sol[eidx]` records the solution.
6. **`remaining_range_for_vertex(v_idx, next_pos)`**: For bundles in `order[next_pos:]`, if the vertex is the `u` end then contribution is `+y ∈ [0,k]`, if the `v` end then `-y ∈ [-k,0]`; accumulate to get the reachable interval `[L, U]`.
7. **`residue_is_feasible(v_idx, next_pos)`**: Let `need = (target_residue[v_idx] - partial_sum[v_idx]) % l`; we need some $t\in [L,U]$ with $t\equiv \text{need}\pmod l$, which is equivalent to some integer $q$ with $L\le \text{need}+ql\le U$; `qmin = ceil_div(L - need, l)`, `qmax = (U - need) // l`; return `qmin <= qmax`.
8. **DFS recursion**: If `pos == len(order)`, check all vertices satisfy `partial_sum[v_idx] % l == target_residue[v_idx]`; otherwise for the current bundle `eidx` try each `y in domains[eidx]`, update `partial_sum`, and if all vertices pass `residue_is_feasible` and recursion succeeds, record `y_sol[eidx]=y` and return; otherwise rollback.
9. **Assembly**: From `y_sol` compute `y_by_pair`, `out_minus_in`, `directions` (first $y$ edges as $u\to v$, remaining $k-y$ as $v\to u$), then verify `out_minus_in[v] % mod2l == beta[v]`.

### 5.5 `is_SZl` and `get_all_infeasible_betas`

`is_SZl`: enumerate all valid $\beta$; if any is infeasible, return `(False, beta)`; otherwise `(True, None)`. `get_all_infeasible_betas`: collect all infeasible $\beta$ into a list and return it.

---

## 6. `generate_nonisomorphic_szl.py`: what it does

This script, for **given $n$, $m$, $l$** (all configurable via command-line arguments):

1. enumerates $n$-vertex, $m$-edge (with multiplicity) multigraphs under constraints: min degree $\ge l-1$, max multiplicity per pair $\le l-2$, edge connectivity $\ge 2$, connected
2. excludes graphs containing known $SZ_l$ subgraphs: 3-vertex (with $\ge 2l-2$ edges), 4-vertex (when $n>4$, by the paper's theorem)
3. removes isomorphic duplicates (keeps representatives)
4. draws overview figures (up to 5×5 graphs per image; multiple images if needed) and saves them in the output folder
5. tests each representative for $SZ_l$ using `SZlSolver(..., l)`; for graphs that are not $SZ_l$, prints a summary of indices first, then **all** infeasible $\beta$ for each

---

## 7. `generate_nonisomorphic_szl.py`: implementation workflow (step-by-step)

### 7.1 Enumeration space and recursive generation

- **`make_pairs(n)`**: returns `[(i,j) for i in 1..n for j in i+1..n]`, i.e. all unordered pairs in lexicographic order.
- **`_enumerate_multiplicities_rec(pairs, remaining, max_per, prefix)`**:
  - If `pairs` is empty, yield `tuple(prefix)` when `remaining==0`, else return;
  - Otherwise for `k` from `min(remaining, max_per)` down to 0, recursively call `_enumerate_multiplicities_rec(pairs[1:], remaining-k, max_per, prefix+[k])`;
  - This enumerates all tuples with `sum=k` and each entry $\le$ max_per.
- **`enumerate_nonisomorphic_graphs`**:
  1. Call `_enumerate_multiplicities_rec(pairs, m, l-2, [])` to get all tuples with total edges $m$ and per-pair $\le l-2$;
  2. For each `counts`: check min degree, connectivity, edge connectivity $\ge 2$, no forbidden 3-vertex subgraph, no forbidden 4-vertex subgraph (when $n>4$);
  3. Deduplicate using `canonical_key(counts)`; `reps[key] = counts`; return `list(reps.values())`.

### 7.2 Subgraph exclusion (step-by-step)

- **`has_forbidden_3vertex_subgraph`**: Use `pair_to_count` to convert `counts` to a `(u,v)->w` dict; triple loop over all 3-vertex sets `(i,j,k)`, compute `m_ij + m_ik + m_jk`; return True if $\ge 2l-2$.
- **`has_forbidden_4vertex_subgraph`**: For each 4-vertex subset `quad`, extract the 6 pair multiplicities and apply SZ_l-simplification (`min(w, l-1)`); call `_is_4vertex_szl_by_theorem(mult_4, l_val)` to decide if that 4-vertex graph is SZ_l; return True if any such subgraph exists.
- **`_is_4vertex_szl_by_theorem`**: Check $\delta\ge l-1$, $e\ge 3l-3$; check (1) $\mu=l-1$, or (2) some 2-vertex cut $\ne 2l-2$, or (3) some vertex degree $\not\equiv l\pmod 2$; return True if any condition holds.

### 7.3 Isomorphism reduction: canonical key (step-by-step)

- **`canonical_key(counts, pairs, n)`**: Convert `counts` to `w[(u,v)] = multiplicity`; enumerate all $n!$ permutations `perm`, for each build inverse map `inv` (new label -> old label); for each pair `(a,b)` in fixed order `order_pairs`, use `inv` to get original pair `(i_old, j_old)` and look up `w[key]`; form the relabeled tuple `tup_t`; keep the lexicographically smallest as `best`; return `best`, used with `reps.setdefault(best, counts)` for deduplication.

### 7.4 Plotting and output (step-by-step)

- **`draw_case(ax, counts, pairs, n, ...)`**: Build simple graph with `nx.Graph` (presence only); use radial layout with center `hub` and others evenly on a circle; `symmetric_rads(m, step)` generates bend radii for $m$ parallel edges (odd includes 0, even does not); for each pair `(u,v)` with `mult` edges, draw each with `FancyArrowPatch` and `arc3,rad=r`; fix `xlim/ylim` for consistent visual scale.
- **Figure output in main**: At most 25 subplots per figure (5×5); if `len(reps) > 25`, split into multiple files (`_overview_1.png`, `_overview_2.png`, ...); for each figure, `fig, axes = plt.subplots(fig_rows, cols, ...)`, call `draw_case` per subplot, then `fig.savefig`, `plt.close`.

### 7.5 Batch $SZ_l$ testing and output (step-by-step)

1. For each representative `counts`, build `nx.MultiGraph` with `build_multigraph(counts, pairs, n)`;
2. `solver = SZlSolver(Gm, l_val)`, `infeasible_betas = solver.get_all_infeasible_betas(verbose=False)`;
3. If infeasible $\beta$ exist, add to `non_sz_cases`;
4. First output `Graphs that are NOT SZ_l: #i1, #i2, ...`;
5. Then for each non-$SZ_l$ graph output: edge multiplicities, degrees, all infeasible $\beta$ as vectors, and total count.

---

## 8. Performance and reproducibility notes

- **Solver complexity**: the number of valid boundaries is roughly $l^{n-1}$ (subject to parity); each DFS is exponential in the worst case. This implementation targets small-to-moderate instances for experimental verification.
- **Determinism**: vertices and bundles are sorted; isomorphism reduction uses the full canonical-key approach, so output is deterministic.
- **Comparison with odd-modulus version**: the odd-modulus scripts can be used when $l$ is odd, with 2 invertible for slightly simpler logic; the general-modulus version handles both odd and even $l$ uniformly.

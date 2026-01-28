from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import importlib.util
import sys
import os

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np


# ---------- 动态加载 SZlSolver ----------

def load_szl_solver_class() -> type:
    here = os.path.dirname(os.path.abspath(__file__))
    solver_path = os.path.join(here, "SZl-identify.py")
    spec = importlib.util.spec_from_file_location("szl_identify_module", solver_path)
    if spec is None or spec.loader is None:
        raise ImportError("无法加载 SZl-identify.py")
    module = importlib.util.module_from_spec(spec)
    # 确保在执行前注册到 sys.modules，便于 dataclasses/typing 解析注解
    sys.modules[spec.name] = module  # type: ignore[arg-type]
    spec.loader.exec_module(module)  # type: ignore[attr-defined]
    if not hasattr(module, "SZlSolver"):
        raise ImportError("SZl-identify.py 中未找到 SZlSolver")
    return module.SZlSolver  # type: ignore[return-value]


# ---------- 组合生成与同构判重 ----------

Pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]


def degree_from_counts(k: Tuple[int, int, int, int, int, int]) -> Dict[int, int]:
    k12, k13, k14, k23, k24, k34 = k
    return {
        1: k12 + k13 + k14,
        2: k12 + k23 + k24,
        3: k13 + k23 + k34,
        4: k14 + k24 + k34,
    }


def is_connected_from_counts(k: Tuple[int, int, int, int, int, int]) -> bool:
    G = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4])
    for (u, v), w in zip(Pairs, k):
        if w > 0:
            G.add_edge(u, v)
    return nx.is_connected(G)


def canonical_key(k: Tuple[int, int, int, int, int, int]) -> Tuple[int, ...]:
    """在所有 4! 顶点置换下，选择字典序最小的上三角 6 元组作为标准形。
    做法：对每个置换 sigma（旧→新），令 inv = sigma^{-1}（新→旧）。
    则新图在 (1,2),(1,3),... 的重数分别为原图 inv(1),inv(2) 等顶点对的重数。
    依此构造 6 元组，取其中字典序最小者。
    """
    labels = [1, 2, 3, 4]
    # 原图的权重字典
    w: Dict[Tuple[int, int], int] = {}
    for (u, v), val in zip(Pairs, k):
        a, b = (u, v) if u < v else (v, u)
        w[(a, b)] = val

    best: Tuple[int, ...] | None = None
    order_pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    for perm in itertools.permutations(labels):
        inv = {perm[i]: labels[i] for i in range(4)}  # 新→旧
        tup = []
        for (a, b) in order_pairs:  # 新标号下的固定顺序
            i_old = inv[a]
            j_old = inv[b]
            key = (i_old, j_old) if i_old < j_old else (j_old, i_old)
            tup.append(w[key])
        tup_t = tuple(tup)
        if best is None or tup_t < best:
            best = tup_t
    assert best is not None
    return best


def enumerate_nonisomorphic_graphs() -> List[Tuple[int, int, int, int, int, int]]:
    reps: Dict[Tuple[int, ...], Tuple[int, int, int, int, int, int]] = {}
    # 每对重数 ∈ [0,3]，总和为 12，且各顶点度 ≥ 4，且连通
    for k12 in range(0, 4):
        for k13 in range(0, 4):
            for k14 in range(0, 4):
                for k23 in range(0, 4):
                    for k24 in range(0, 4):
                        for k34 in range(0, 4):
                            k = (k12, k13, k14, k23, k24, k34)
                            if sum(k) != 12:
                                continue
                            deg = degree_from_counts(k)
                            if min(deg.values()) < 4:
                                continue
                            if not is_connected_from_counts(k):
                                continue
                            key = canonical_key(k)
                            reps.setdefault(key, k)
    return list(reps.values())


# ---------- 绘图（小图、每行 5 个） ----------

def symmetric_rads(m: int, step: float) -> List[float]:
    if m <= 1:
        return [0.0]
    start = -step * (m - 1) / 2.0
    return [start + i * step for i in range(m)]


def draw_case(ax, counts: Tuple[int, int, int, int, int, int], *, hub: int = 1, top: int | None = 2,
              radius: float = 0.8, node_size: int = 280, curve_step: float = 0.18,
              node_color: str = "#66D1B3", edge_color: str = "#5AA9E6") -> None:
    # 节点与径向坐标
    H = nx.Graph()
    H.add_nodes_from([1, 2, 3, 4])
    for (u, v), w in zip(Pairs, counts):
        if w > 0:
            H.add_edge(u, v)
    pos = {hub: (0.0, 0.0)}
    others = [v for v in [1, 2, 3, 4] if v != hub]
    # 旋转使 top 在正上（若给定且存在）
    if top in others:
        idx = others.index(top)  # type: ignore[arg-type]
        others = others[idx:] + others[:idx]
    import math
    for i, v in enumerate(others):
        theta = math.pi / 2 + 2 * math.pi * i / len(others)
        pos[v] = (radius * math.cos(theta), radius * math.sin(theta))

    # 节点与标签
    nx.draw_networkx_nodes(H, pos, node_color=node_color, node_size=node_size, linewidths=1.0, edgecolors="#2C3E50", ax=ax)
    nx.draw_networkx_labels(H, pos, font_color="#000000", font_size=9, ax=ax)
    ax.set_aspect('equal', adjustable='datalim')

    # 多重边为弯曲曲线
    for (u, v), m in zip(Pairs, counts):
        if m == 0:
            continue
        rads = symmetric_rads(m, curve_step)
        for rad in rads:
            # 奇数含直线；偶数全弯曲
            r = rad
            if m % 2 == 0 and abs(r) < 1e-6:
                r = 0.12
            patch = FancyArrowPatch(pos[u], pos[v], connectionstyle=f"arc3,rad={r}", arrowstyle='-',
                                    color=edge_color, linewidth=1.4, alpha=0.95, zorder=1)
            ax.add_patch(patch)

    # 边距与轴
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    pad = 0.22
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.axis('off')


# ---------- 构建 MultiGraph ----------

def build_multigraph(counts: Tuple[int, int, int, int, int, int]) -> nx.MultiGraph:
    Gm = nx.MultiGraph()
    Gm.add_nodes_from([1, 2, 3, 4])
    for (u, v), w in zip(Pairs, counts):
        for _ in range(w):
            Gm.add_edge(u, v)
    return Gm


# ---------- 主流程 ----------

def main():
    l_val = 5
    reps = enumerate_nonisomorphic_graphs()
    print(f"满足条件的非同构图数量: {len(reps)}")

    # 绘图：每行 5 个
    cols = 5
    n = len(reps)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.2, rows * 2.2))
    # 扁平化 axes（稳妥处理 ndarray/list/单轴）
    if isinstance(axes, np.ndarray):
        axes_flat = list(axes.flat)
    elif isinstance(axes, (list, tuple)):
        axes_flat: List = []
        for row in axes:
            if isinstance(row, (list, tuple, np.ndarray)):
                axes_flat.extend(list(np.array(row).flat))
            else:
                axes_flat.append(row)
    else:
        axes_flat = [axes]

    for i, counts in enumerate(reps):
        ax = axes_flat[i]
        draw_case(ax, counts, hub=1, top=2, radius=0.8)
        ax.set_title(f"#{i+1}", fontsize=9)

    # 清理多余轴
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis('off')

    plt.tight_layout()
    plt.show()

    # SZ_l 检测
    SZlSolver = load_szl_solver_class()
    non_sz_cases: List[Tuple[int, Tuple[int, int, int, int, int, int]]] = []
    for i, counts in enumerate(reps):
        Gm = build_multigraph(counts)
        solver = SZlSolver(Gm, l_val)
        ok, witness = solver.is_SZl(verbose=False)
        if not ok:
            non_sz_cases.append((i, counts))
            print("\n图 #{} 不为 SZ_{}:".format(i + 1, l_val))
            print("- 边重数 (12 总): {}".format(dict(zip(Pairs, counts))))
            deg = degree_from_counts(counts)
            print("- 顶点度: {}".format(deg))
            if witness is not None:
                print("- 一个不可实现的 beta（向量表示）：{}".format(solver.format_beta(witness)))

    if not non_sz_cases:
        print("\n所有图均为 SZ_{}。".format(l_val))


if __name__ == "__main__":
    # 该脚本生成所有满足条件的非同构多重图、展示小图网格，并对每个图检测是否为 SZ_5。
    main()



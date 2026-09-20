"""
The two-way split heatmap.

This is the figure ``ComplexHeatmap`` draws when you hand it ``row_split`` and
``column_split`` from k-means and leave ``cluster_rows`` and ``cluster_columns`` on:
k-means decides the **blocks**, hierarchical clustering decides the **order inside
each block**, and every block gets its own small dendrogram. k-means is a partition,
not a tree, so the blocks carry no nesting and no branch lengths; the trees you see
are local to a block and say nothing about how the blocks relate to one another.

Both axes are always clustered. Rows are samples, columns are the measured objects,
and the default pairing is minkowski with ward.D2 on both.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

MISSING = {"NA", "nan", "NaN", "None", "", "<NA>"}


def _block_order(M: pd.DataFrame, axis: str, assign: pd.Series,
                 method_dist: str, method_hclust: str):
    """Order within each k-means block by hierarchical clustering.

    Returns ``(labels_in_order, blocks)`` where each block is
    ``(cluster_id, [labels], linkage_or_None)``.
    """
    from pvclust_py.distance import distance
    from pvclust_py.hclust import linkage
    from scipy.cluster.hierarchy import dendrogram as sd

    order, blocks = [], []
    for cid in sorted(assign.unique()):
        members = [l for l in assign.index[assign == cid]]
        sub = M.loc[members] if axis == "rows" else M[members]
        Z = None
        if len(members) > 2:
            # distance() clusters the COLUMNS of what it is given, so a row block
            # has to be transposed into that orientation first.
            A = sub.to_numpy(float).T if axis == "rows" else sub.to_numpy(float)
            D = distance(A, method_dist)
            if np.isfinite(D).all():
                Z = linkage(D, method_hclust)
                members = [members[i] for i in sd(Z, no_plot=True)["leaves"]]
        order += members
        blocks.append((cid, members, Z))
    return order, blocks


def _strip(ax, ann: pd.DataFrame, horizontal: bool):
    """Annotation tiles. Each level of each variable gets its own colour, and the
    mapping is returned so the caller can draw a key -- sharing one colormap across
    variables makes the same colour mean 'female' in one strip and 'remission' in the
    next, which is worse than no colour at all."""
    import matplotlib.pyplot as plt

    palette = (list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
               + list(plt.get_cmap("tab20c").colors))
    grey = (0.88, 0.88, 0.88)
    rows, key, k = [], [], 0
    for col in ann.columns:
        s = ann[col].astype(str)
        levels = sorted({v for v in s.unique() if v not in MISSING})
        colours = {}
        for lv in levels:
            colours[lv] = palette[k % len(palette)]
            k += 1
        key.append((col, [(lv, colours[lv]) for lv in levels]))
        rows.append([colours.get(v, grey) for v in s])
    A = np.array(rows, dtype=float)
    ax.imshow(A if horizontal else A.transpose(1, 0, 2), aspect="auto",
              interpolation="nearest")
    if horizontal:
        ax.set_yticks(range(len(ann.columns)))
        ax.set_yticklabels(ann.columns, fontsize=6)
        ax.set_xticks([])
    else:
        ax.set_xticks(range(len(ann.columns)))
        ax.set_xticklabels(ann.columns, rotation=90, fontsize=6)
        ax.set_yticks([])
    return key


def _block_trees(ax, blocks, horizontal: bool):
    """One small dendrogram per block, laid over that block's span."""
    from scipy.cluster.hierarchy import dendrogram as sd

    ax.set_xticks([]); ax.set_yticks([])
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)

    total = sum(len(m) for _cid, m, _Z in blocks)
    heights = [float(Z[:, 2].max()) for _c, _m, Z in blocks if Z is not None and Z.size]
    hmax = max(heights) if heights else 1.0
    start = 0
    for _cid, members, Z in blocks:
        n = len(members)
        if Z is not None and Z.size:
            dd = sd(Z, no_plot=True)
            for xs, ys in zip(dd["icoord"], dd["dcoord"]):
                # scipy lays leaves at 5, 15, 25...; map onto this block's own span
                px = [start + (x - 5.0) / 10.0 + 0.5 for x in xs]
                py = [y / hmax for y in ys]
                if horizontal:
                    ax.plot(px, py, color="#555555", lw=0.9)
                else:
                    ax.plot([1.0 - v for v in py], px, color="#555555", lw=0.9)
        start += n
    if horizontal:
        ax.set_xlim(0, total); ax.set_ylim(0, 1.05)
    else:
        ax.set_ylim(total, 0); ax.set_xlim(0, 1.05)


def split_heatmap(matrix, base: str, *, row_assign: pd.Series, col_assign: pd.Series,
                  row_annotations: Optional[pd.DataFrame] = None,
                  method_dist: str = "minkowski", method_hclust: str = "ward.D2",
                  z_score: Optional[str] = "columns", cmap: str = "RdBu_r",
                  title: Optional[str] = None, max_labels: int = 60) -> None:
    """Heatmap split into k-means blocks on both axes -> png, svg, pdf and html.

    Args:
        row_assign: cluster id per row, from :func:`kmeans_py.core.kmeans`.
        col_assign: cluster id per column.
        row_annotations: DataFrame indexed like the rows, one tile row per column.
        method_dist, method_hclust: how the order WITHIN each block is decided.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    M = matrix if isinstance(matrix, pd.DataFrame) else pd.DataFrame(matrix)
    row_assign = row_assign.reindex(M.index)
    col_assign = col_assign.reindex(M.columns)

    row_order, row_blocks = _block_order(M, "rows", row_assign, method_dist, method_hclust)
    col_order, col_blocks = _block_order(M, "columns", col_assign, method_dist, method_hclust)
    M = M.loc[row_order, col_order]

    if z_score == "columns":
        M = (M - M.mean()) / M.std().replace(0, np.nan)
    elif z_score == "rows":
        M = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
    M = M.fillna(0.0)

    n_ann = 0 if row_annotations is None else row_annotations.shape[1]
    w_cap = 40 if M.shape[1] <= max_labels else 26
    h_cap = 40 if M.shape[0] <= max_labels else 24
    fig = plt.figure(figsize=(max(9, min(w_cap, 0.16 * M.shape[1] + 5)),
                              max(7, min(h_cap, 0.14 * M.shape[0] + 4))),
                     constrained_layout=True)
    gs = GridSpec(3, 3 + (1 if n_ann else 0), figure=fig,
                  height_ratios=[1.1, 6, 0.35],
                  width_ratios=[1.0] + ([0.12 * n_ann] if n_ann else []) + [6, 0.25])
    c_main = 1 + (1 if n_ann else 0)

    _block_trees(fig.add_subplot(gs[0, c_main]), col_blocks, horizontal=True)
    _block_trees(fig.add_subplot(gs[1, 0]), row_blocks, horizontal=False)

    key = []
    if n_ann:
        key = _strip(fig.add_subplot(gs[1, 1]),
                     row_annotations.reindex(M.index), horizontal=False)

    ax = fig.add_subplot(gs[1, c_main])
    lim = float(np.nanpercentile(np.abs(M.to_numpy()), 99)) or 1.0
    im = ax.imshow(M.to_numpy(), aspect="auto", cmap=cmap, vmin=-lim, vmax=lim,
                   interpolation="nearest")

    # white rules between blocks -- what makes the split readable
    x = 0
    for _cid, members, _Z in col_blocks[:-1]:
        x += len(members)
        ax.axvline(x - 0.5, color="white", lw=2.5)
    y = 0
    for _cid, members, _Z in row_blocks[:-1]:
        y += len(members)
        ax.axhline(y - 0.5, color="white", lw=2.5)

    if M.shape[1] <= max_labels:
        # Names go at the TOP, directly under the column dendrogram, so a branch and
        # the name it belongs to can be read together. At the bottom of a tall figure
        # the tree and its labels are hundreds of rows apart and the pairing is lost.
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels(col_order, rotation=90, fontsize=6)
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("bottom")
    else:
        ax.set_xticks([])
    if M.shape[0] <= max_labels:
        ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(row_order, fontsize=6)
    else:
        ax.set_yticks([])
    ax.set_xlabel(f"{M.shape[1]} columns in {len(col_blocks)} k-means blocks"
                  + ("" if M.shape[1] <= max_labels else " (labels suppressed)"))
    ax.set_ylabel(f"{M.shape[0]} rows in {len(row_blocks)} k-means blocks"
                  + ("" if M.shape[0] <= max_labels else " (labels suppressed)"))

    cax = fig.add_subplot(gs[1, c_main + 1])
    fig.colorbar(im, cax=cax, label="z-score" if z_score else "value")

    if key:
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        handles, labels = [], []
        for col, levels in key:
            if not levels:
                continue
            handles.append(Line2D([], [], linestyle="none"))
            labels.append(f"$\\bf{{{col.replace('_', chr(92) + '_')}}}$")
            if len(levels) > 12:
                handles.append(Line2D([], [], linestyle="none"))
                labels.append(f"  {len(levels)} levels")
                continue
            for lv, colour in levels:
                handles.append(Patch(facecolor=colour, edgecolor="none"))
                labels.append(f"  {lv}")
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
                   frameon=False, fontsize=7, handlelength=1.2, borderaxespad=0.0)

    if title:
        fig.suptitle(title, fontsize=12, fontweight="bold")
    # The names as text too: a long list is easier to search than to squint at.
    pd.DataFrame([{"block": cid, "position": i, "object": m}
                  for cid, members, _Z in col_blocks
                  for i, m in enumerate(members)]).to_csv(base + "_column_blocks.csv",
                                                          index=False)
    pd.DataFrame([{"block": cid, "position": i, "object": m}
                  for cid, members, _Z in row_blocks
                  for i, m in enumerate(members)]).to_csv(base + "_row_blocks.csv",
                                                          index=False)

    fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(base + ".svg", bbox_inches="tight")
    # PDF as well: vector, so it prints at any paper size without pixelating, and it
    # is the format a printer or a journal will actually accept.
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)

    _interactive(M, base, title, row_order, col_order, row_blocks, col_blocks,
                 row_annotations)


def _interactive(M, base, title, row_order, col_order, row_blocks, col_blocks, ann):
    """Self-contained plotly version, block membership on hover."""
    try:
        import plotly.graph_objects as go
    except Exception:
        return
    rblock = {m: str(cid) for cid, members, _Z in row_blocks for m in members}
    cblock = {m: str(cid) for cid, members, _Z in col_blocks for m in members}
    text = np.empty(M.shape, dtype=object)
    for i, r in enumerate(row_order):
        extra = ""
        if ann is not None and r in ann.index:
            extra = "<br>" + "<br>".join(f"{c}: {ann.loc[r, c]}" for c in ann.columns)
        for j, c in enumerate(col_order):
            text[i, j] = (f"row {r} (block {rblock.get(r, '?')})<br>"
                          f"col {c} (block {cblock.get(c, '?')})<br>"
                          f"z = {M.iat[i, j]:.2f}{extra}")
    fig = go.Figure(go.Heatmap(z=M.to_numpy(), x=list(col_order), y=list(row_order),
                               colorscale="RdBu_r", reversescale=False,
                               text=text, hoverinfo="text"))
    fig.update_layout(title_text=title or base,
                      width=min(1500, 60 + 14 * M.shape[1]),
                      height=min(1400, 120 + 12 * M.shape[0]),
                      yaxis=dict(autorange="reversed"))
    fig.write_html(base + ".html")

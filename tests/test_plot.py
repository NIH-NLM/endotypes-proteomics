"""The two-way split heatmap.

k-means decides the blocks; hierarchical clustering orders within each block and draws
its own small dendrogram. Both axes, always.
"""
import numpy as np
import pandas as pd
import pytest

from kmeans_py import kmeans
from kmeans_py.plot import _block_order, split_heatmap


@pytest.fixture
def cohort():
    rng = np.random.default_rng(0)
    blocks = []
    for g, shift in enumerate([0.0, 2.5, -2.5]):
        blocks.append(pd.DataFrame(rng.normal(shift, 0.6, size=(20, 12)),
                                   index=[f"g{g}_s{i}" for i in range(20)]))
    X = pd.concat(blocks)
    X.columns = [f"p{j}" for j in range(12)]
    meta = pd.DataFrame({"Group": ["HV" if i % 3 else "SLE" for i in range(len(X))],
                         "Sex": ["F" if i % 4 else "M" for i in range(len(X))]},
                        index=X.index)
    return X, meta


def test_block_order_keeps_every_label_once(cohort):
    X, _ = cohort
    a = kmeans(X, 3, cluster="rows").assignment
    order, blocks = _block_order(X, "rows", a, "minkowski", "ward.D2")
    assert sorted(order) == sorted(X.index)
    assert len(blocks) == 3
    assert sorted(m for _c, ms, _Z in blocks for m in ms) == sorted(X.index)


def test_block_order_does_not_interleave_blocks(cohort):
    """Members of one block must be contiguous, or the white rules land wrong."""
    X, _ = cohort
    a = kmeans(X, 3, cluster="rows").assignment
    order, blocks = _block_order(X, "rows", a, "minkowski", "ward.D2")
    pos = {lab: i for i, lab in enumerate(order)}
    for _cid, members, _Z in blocks:
        idx = sorted(pos[m] for m in members)
        assert idx == list(range(idx[0], idx[0] + len(idx))), "block is not contiguous"


def test_each_block_gets_its_own_tree(cohort):
    X, _ = cohort
    a = kmeans(X, 3, cluster="rows").assignment
    _order, blocks = _block_order(X, "rows", a, "minkowski", "ward.D2")
    for _cid, members, Z in blocks:
        if len(members) > 2:
            assert Z is not None and Z.shape == (len(members) - 1, 4)


def test_writes_all_three_formats(cohort, tmp_path):
    X, meta = cohort
    rows = kmeans(X, 3, cluster="rows")
    cols = kmeans(X, 2, cluster="columns")
    base = str(tmp_path / "hm")
    split_heatmap(X, base, row_assign=rows.assignment, col_assign=cols.assignment,
                  row_annotations=meta)
    for ext in (".png", ".svg", ".html"):
        assert (tmp_path / f"hm{ext}").exists(), ext
        assert (tmp_path / f"hm{ext}").stat().st_size > 0


def test_both_axes_are_split(cohort, tmp_path):
    """The figure must reflect BOTH partitions, not just the rows."""
    X, meta = cohort
    rows = kmeans(X, 3, cluster="rows")
    cols = kmeans(X, 2, cluster="columns")
    _o, rb = _block_order(X, "rows", rows.assignment, "minkowski", "ward.D2")
    _o, cb = _block_order(X, "columns", cols.assignment, "minkowski", "ward.D2")
    assert len(rb) == 3 and len(cb) == 2

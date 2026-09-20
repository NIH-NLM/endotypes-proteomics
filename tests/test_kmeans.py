"""The partition, and the description of it."""
import numpy as np
import pandas as pd
import pytest

from kmeans_py import associations, choose_k, composition, distinguishing, kmeans


@pytest.fixture
def cohort():
    """Three groups of samples with genuinely different profiles, plus metadata."""
    rng = np.random.default_rng(0)
    blocks, meta = [], []
    for g, (shift, label) in enumerate([(0.0, "HV"), (2.5, "SLE"), (-2.5, "SLE")]):
        n = 30
        blocks.append(pd.DataFrame(rng.normal(shift, 0.6, size=(n, 12)),
                                   index=[f"g{g}_s{i}" for i in range(n)]))
        meta += [{"sample": f"g{g}_s{i}", "Group": label,
                  "Sex": "F" if i % 4 else "M", "SLEDAI": float(g * 3 + i % 4)}
                 for i in range(n)]
    X = pd.concat(blocks)
    X.columns = [f"p{j}" for j in range(12)]
    return X, pd.DataFrame(meta).set_index("sample")


def test_partitions_every_sample_once(cohort):
    X, _ = cohort
    res = kmeans(X, k=3)
    assert len(res.assignment) == len(X)
    assert sorted(res.assignment.index) == sorted(X.index)
    assert res.sizes.sum() == len(X)


def test_recovers_planted_groups(cohort):
    X, _ = cohort
    res = kmeans(X, k=3)
    for members in res.members.values():
        assert len({m.split("_")[0] for m in members}) == 1


def test_clusters_the_requested_axis(cohort):
    X, _ = cohort
    assert set(kmeans(X, k=3, cluster="rows").assignment.index) == set(X.index)
    assert set(kmeans(X, k=3, cluster="columns").assignment.index) == set(X.columns)


def test_centroid_tree_has_one_leaf_per_cluster(cohort):
    X, _ = cohort
    res = kmeans(X, k=4)
    assert res.linkage is not None
    assert res.linkage.shape == (3, 4)          # k - 1 merges


def test_silhouette_reports_each_cluster_and_overall(cohort):
    X, _ = cohort
    res = kmeans(X, k=3)
    assert "overall" in res.silhouette.index
    assert res.silhouette["overall"] > 0.5      # planted groups are well separated


def test_k_must_be_reachable(cohort):
    X, _ = cohort
    with pytest.raises(ValueError, match="k must be between"):
        kmeans(X, k=len(X) + 1)
    with pytest.raises(ValueError, match="cluster must be"):
        kmeans(X, k=2, cluster="diagonal")


def test_choose_k_peaks_at_the_planted_number(cohort):
    X, _ = cohort
    out = choose_k(X, range(2, 7))
    assert int(out.loc[out["silhouette"].idxmax(), "k"]) == 3


def test_composition_counts_each_cluster(cohort):
    X, meta = cohort
    res = kmeans(X, k=3)
    comp = composition(res.assignment, meta, numeric=["SLEDAI"])
    group = comp[comp["variable"] == "Group"]
    for cid, sub in group.groupby("cluster"):
        assert sub["n"].sum() == res.sizes[int(cid)]
        assert abs(sub["pct"].sum() - 100.0) < 1e-6
    sledai = comp[comp["variable"] == "SLEDAI"]
    assert sledai["median"].notna().all()
    assert sledai["level"].eq("").all()         # numeric variables have no levels


def test_associations_find_the_planted_signal(cohort):
    X, meta = cohort
    res = kmeans(X, k=3)
    assoc = associations(res.assignment, meta, numeric=["SLEDAI"]).set_index("variable")
    assert assoc.loc["Group", "p_value"] < 0.05
    assert assoc.loc["Sex", "p_value"] > 0.05   # planted independent of the groups


def test_distinguishing_is_signed_and_ranked(cohort):
    X, _ = cohort
    res = kmeans(X, k=3)
    d = distinguishing(X, res.assignment, top=4)
    assert len(d) == 3 * 4
    for _cid, sub in d.groupby("cluster"):
        mags = sub["std_mean_difference"].abs().to_numpy()
        assert (np.diff(mags) <= 1e-9).all()    # ranked by magnitude
    assert set(d["direction"]) <= {"higher", "lower"}

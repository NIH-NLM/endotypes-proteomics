"""The federated partition.

k-medoids, because the pooled coordinates must never exist anywhere: Lloyd's algorithm
averages points into centroids and so needs them, while medoids need only distances,
and pvclust-py recovers the pooled distance matrix exactly from sufficient statistics.
"""
import numpy as np
import pandas as pd
import pytest

from pvclust_py.distance import distance

from kmeans_py import federated_kmeans, kmeans


@pytest.fixture
def blocks():
    """Four tight blocks of five objects -- a partition any method should find."""
    rng = np.random.default_rng(0)
    cols = {}
    for b in range(4):
        driver = rng.normal(size=80)
        for j in range(5):
            cols[f"b{b}_{j}"] = driver + rng.normal(scale=0.2, size=80)
    return pd.DataFrame(cols)


def test_partition_covers_every_object_once(blocks):
    D = distance(blocks.to_numpy(float), "correlation")
    cat = federated_kmeans(D, list(blocks.columns), k=4)
    assert len(cat) == 4
    assert sorted(m for e in cat for m in e["members"]) == sorted(blocks.columns)


def test_partition_recovers_the_blocks(blocks):
    D = distance(blocks.to_numpy(float), "correlation")
    for e in federated_kmeans(D, list(blocks.columns), k=4):
        assert len({m.split("_")[0] for m in e["members"]}) == 1, e["members"]


def test_pooled_partition_matches_the_local_one(blocks):
    """The point of federating: distances alone give the answer coordinates would."""
    D = distance(blocks.to_numpy(float), "correlation")
    federated = {frozenset(e["members"]) for e in federated_kmeans(D, list(blocks.columns), k=4)}
    local = {frozenset(m) for m in kmeans(blocks, k=4, cluster="columns").members.values()}
    assert federated == local


def test_k_must_fit_the_object_count(blocks):
    D = distance(blocks.to_numpy(float), "correlation")
    with pytest.raises(ValueError, match="k must be between"):
        federated_kmeans(D, list(blocks.columns), k=len(blocks.columns) + 1)


def test_distance_matrix_must_match_the_labels(blocks):
    D = distance(blocks.to_numpy(float), "correlation")
    with pytest.raises(ValueError, match="distance matrix is"):
        federated_kmeans(D, list(blocks.columns)[:-1], k=3)

"""
Federating a flat partition.

``pvclust-py`` pools sufficient statistics into the exact pooled distance matrix. This
module takes it from there: given distances and nothing else, it finds the partition.

This is the part of the k-means story that federates cleanly, and it holds regardless
of what you think about *p*-values on clusters. Lloyd's algorithm averages points into
centroids, so it needs coordinates -- and the pooled coordinates would be every
project's samples side by side, which is the one thing that must not leave a project.
The objective k-means minimises depends only on pairwise distances (Huygens' theorem),
and the pooled distance matrix is recoverable exactly from sufficient statistics. So
the partition can be found centrally without the coordinates existing anywhere.
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np


def federated_kmeans(D, labels: Sequence[str], k: int, *, max_iter: int = 300,
                     n_init: int = 10, seed: int = 42) -> List[Dict]:
    """A federated k-means-style partition, from the pooled distance matrix alone.

    Why this works. Ordinary k-means (Lloyd's) needs coordinates, to average points
    into centroids -- and the pooled coordinate matrix would be everyone's samples
    side by side, which is exactly what must not leave a project. But the objective
    k-means minimises, the within-cluster sum of squares, depends only on PAIRWISE
    DISTANCES (Huygens' theorem). Since :func:`pvclust_py.aggregate.federated_tree` already recovers the
    pooled distance matrix exactly from sufficient statistics, the partition can be
    found centrally without any coordinates existing anywhere.

    This is k-medoids (PAM-style): the cluster representative is an actual object
    rather than an average, so nothing needs to be averaged. On euclidean-type
    distances it optimises the same objective as k-means; the medoid constraint makes
    it slightly more restrictive and, as a rule, more robust to outliers.

    Args:
        D: pooled ``p x p`` distance matrix from :func:`pvclust_py.aggregate.federated_tree`.
        labels: object names, in the same order as ``D``.
        k: number of clusters.

    Returns:
        Edge records in the same shape :func:`~pvclust_py.hclust.edge_table` produces,
        so they can be used as a catalogue for ``count-edges`` exactly like tree edges.
    """
    from pvclust_py.hclust import edge_id

    D = np.asarray(D, dtype=float)
    p = len(labels)
    if D.shape != (p, p):
        raise ValueError(f"distance matrix is {D.shape}, expected ({p}, {p})")
    if not 2 <= k <= p:
        raise ValueError(f"k must be between 2 and {p}, got {k}")

    rng = np.random.default_rng(seed)
    best_cost, best_assign = np.inf, None

    for _ in range(n_init):
        # k-means++ style seeding, on distances rather than coordinates
        medoids = [int(rng.integers(p))]
        for _ in range(k - 1):
            d2 = D[medoids].min(axis=0) ** 2
            total = d2.sum()
            medoids.append(int(rng.choice(p, p=d2 / total)) if total > 0
                           else int(rng.integers(p)))
        medoids = np.array(medoids)

        for _ in range(max_iter):
            assign = np.argmin(D[medoids], axis=0)
            moved = False
            for c in range(k):
                members = np.where(assign == c)[0]
                if len(members) == 0:
                    continue
                # the medoid is the member minimising total distance to the rest
                within = D[np.ix_(members, members)].sum(axis=1)
                candidate = members[int(np.argmin(within))]
                if candidate != medoids[c]:
                    medoids[c], moved = candidate, True
            if not moved:
                break

        cost = float(D[medoids, np.arange(p)].sum()) if False else \
            float(D[medoids].min(axis=0).sum())
        if cost < best_cost:
            best_cost, best_assign = cost, np.argmin(D[medoids], axis=0)

    labels = list(labels)
    out = []
    for c in range(k):
        members = sorted(labels[i] for i in np.where(best_assign == c)[0])
        if not members:
            continue
        out.append({"edge_id": edge_id(members), "members": members,
                    "n_members": len(members), "height": float("nan"),
                    "merge_order": len(out) + 1})
    return out

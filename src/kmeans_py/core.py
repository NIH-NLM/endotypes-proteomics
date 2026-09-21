"""
k-means, plainly.

No *p*-values. The multiscale bootstrap that ``pvclust-py`` applies to a dendrogram
does not transfer to a flat partition: its accuracy comes from an expansion around a
smooth region boundary, and for hierarchical clustering the edges are nested, so small
clusters recur across resamples and the bootstrap probability stays in a useful range.
A k-means cluster has no such nesting. The event "this exact member set is one of the
k clusters" has a probability that vanishes as the object count grows -- measured on
lung expression at k=3, the largest bootstrap probability was 0.825 over 20 objects,
0.202 over 50, and 0.000 over 150. With zero at every scale there is no curve to fit,
so the AU value is not small, it is undefined.

What this package reports instead is what k-means can honestly support: the partition,
how separated it is (silhouette), how each cluster differs from the rest, and how the
clusters line up against the demographics you already hold. For a published
cluster-wise stability measure, see Hennig (2007) and ``fpc::clusterboot``, which
resamples and scores each cluster by its best Jaccard match.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd


@dataclass
class KMeansResult:
    """One k-means partition.

    Attributes:
        assignment: cluster id per clustered entity, indexed by its label.
        centroids: ``k x m`` profile of each cluster, over the *other* axis.
        sizes: members per cluster.
        silhouette: mean silhouette per cluster, and ``overall``. Roughly, how much
            closer a member sits to its own cluster than to the next nearest one.
            Near 1 is well separated, near 0 is on a boundary, negative suggests the
            member belongs elsewhere.
        inertia: within-cluster sum of squares, the quantity k-means minimises. Useful
            only for comparing runs at the same k on the same data.
        linkage: scipy linkage over the CENTROIDS, or None when k < 2. Its leaves are
            the clusters, not the original entities, so it shows how the clusters sit
            relative to one another. k-means itself is flat; this is a picture drawn
            on top of it, and its merge heights carry no inferential claim.
        cluster: which axis was clustered, ``rows`` or ``columns``.
    """
    assignment: pd.Series
    centroids: pd.DataFrame
    sizes: pd.Series
    silhouette: pd.Series
    inertia: float
    linkage: Optional[np.ndarray]
    cluster: str
    k: int

    @property
    def members(self) -> Dict[str, List[str]]:
        """Cluster id -> the labels in it."""
        return {str(c): list(self.assignment.index[self.assignment == c])
                for c in sorted(self.assignment.unique())}

    def table(self) -> pd.DataFrame:
        """One row per cluster: size, silhouette, members."""
        return pd.DataFrame({
            "cluster": [str(c) for c in sorted(self.assignment.unique())],
            "n_members": [self.sizes[c] for c in sorted(self.assignment.unique())],
            "silhouette": [self.silhouette.get(c, np.nan)
                           for c in sorted(self.assignment.unique())],
            "members": [";".join(self.members[str(c)])
                        for c in sorted(self.assignment.unique())],
        })


def kmeans(X, k: int = 3, *, cluster: str = "rows", seed: int = 42,
           n_init: int = 10, centroid_tree: bool = True) -> KMeansResult:
    """Partition one axis of ``X`` into ``k`` clusters.

    Args:
        X: DataFrame, rows x columns.
        cluster: ``rows`` (default) clusters the samples, which is what you want when
            the clusters are to be described by demographics. ``columns`` clusters the
            measured objects -- proteins, genes -- instead.
        n_init: restarts. Do not lower it; k-means is sensitive to initialisation and
            a single start makes the answer depend on the seed more than on the data.
        centroid_tree: also build a dendrogram over the cluster centroids.

    Returns:
        :class:`KMeansResult`.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_samples

    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X)
    if cluster not in ("rows", "columns"):
        raise ValueError(f"cluster must be rows or columns, got {cluster!r}")

    # sklearn clusters rows, so clustering the columns means transposing first.
    A = X if cluster == "rows" else X.T
    labels = list(A.index)
    M = np.nan_to_num(A.to_numpy(dtype=float))
    n = len(labels)
    if not 2 <= k <= n:
        raise ValueError(f"k must be between 2 and the number of clustered items "
                         f"({n}), got {k}")

    km = KMeans(n_clusters=k, n_init=n_init, random_state=seed).fit(M)
    assignment = pd.Series(km.labels_, index=pd.Index(labels, name=A.index.name),
                           name="cluster")

    sil = pd.Series(silhouette_samples(M, km.labels_), index=labels)
    per_cluster = sil.groupby(assignment.to_numpy()).mean()
    per_cluster["overall"] = float(sil.mean())

    centroids = pd.DataFrame(km.cluster_centers_, columns=A.columns,
                             index=[f"cluster{c}" for c in range(k)])

    Z = None
    if centroid_tree and k >= 2:
        from pvclust_py.distance import distance
        from pvclust_py.hclust import linkage as _linkage
        Z = _linkage(distance(centroids.to_numpy(float).T,
                              "correlation" if k > 2 else "euclidean"), "average")

    return KMeansResult(assignment=assignment, centroids=centroids,
                        sizes=assignment.value_counts().sort_index(),
                        silhouette=per_cluster, inertia=float(km.inertia_),
                        linkage=Z, cluster=cluster, k=k)


def choose_k(X, candidates: Sequence[int] = range(2, 11), *, cluster: str = "rows",
             seed: int = 42, n_init: int = 10) -> pd.DataFrame:
    """Mean silhouette and inertia at each k, so the choice is made in the open.

    There is no *p*-value here and no k is "significant". Silhouette peaks at the k
    whose clusters are best separated, and inertia always falls with k, so read the
    two together and say which you used.
    """
    rows = []
    for k in candidates:
        try:
            res = kmeans(X, int(k), cluster=cluster, seed=seed, n_init=n_init,
                         centroid_tree=False)
        except ValueError:
            continue
        rows.append({"k": int(k), "silhouette": float(res.silhouette["overall"]),
                     "inertia": res.inertia,
                     "smallest_cluster": int(res.sizes.min())})
    return pd.DataFrame(rows)

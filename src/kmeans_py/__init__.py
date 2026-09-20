"""
kmeans-py -- k-means clustering, described by the demographics you already hold.

A flat partition and an honest account of it: how separated the clusters are, what
distinguishes each one, and how they line up against disease group, sex, age band or
disease activity. Optionally a dendrogram over the cluster centroids, which shows how
the clusters relate without pretending the partition was hierarchical.

No AU *p*-values. See :mod:`kmeans_py.core` for why the multiscale bootstrap does not
transfer to a flat partition, and Hennig (2007) for the published cluster-wise
stability measure if that is what you need.

Companion to `pvclust-py <https://github.com/NIH-NLM/pvclust-py>`_, which does the
hierarchical case with AU p-values and supplies the shared input handling, batch
correction, distance layer and plotting.
"""
from .aggregate import federated_kmeans
from .core import KMeansResult, choose_k, kmeans
from .profile import associations, composition, distinguishing

__all__ = ["kmeans", "choose_k", "KMeansResult", "federated_kmeans",
           "composition", "associations", "distinguishing"]
__version__ = "0.1.0"

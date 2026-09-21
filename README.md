# kmeans-py

k-means clustering, described by the demographics you already hold.

A flat partition and an honest account of it: how separated the clusters are, what
distinguishes each one, and how they line up against disease group, sex, age band or
disease activity. Optionally a dendrogram over the cluster centroids, which shows how
the clusters relate without pretending the partition was hierarchical.

Companion to [pvclust-py](https://github.com/NIH-NLM/pvclust-py), which does the
hierarchical case with AU *p*-values. That package is a dependency, not a copy: the
input handling, batch correction, distance layer, association statistics and plotting
all live there and are called from here.

## No p-values on the clusters, and why

`pvclust-py` attaches an approximately-unbiased *p*-value to every branch of a
dendrogram. That does not transfer to a flat partition, and we should not pretend it
does.

The multiscale bootstrap gets its accuracy from an expansion around a smooth region
boundary. In hierarchical clustering the region is "this edge appears in the tree", and
edges are nested, so small clusters recur across resamples and the bootstrap
probability stays in a range where the expansion means something. A k-means cluster has
no such nesting. The event is "this exact member set is one of the k clusters", and its
probability vanishes as the object count grows. Measured on lung expression at k=3:

| objects | largest bootstrap probability |
|---|---|
| 20 | 0.825 |
| 50 | 0.202 |
| 150 | 0.000 |

At zero across every scale there is no curve to fit, so the resulting value is not
small, it is undefined. The published route for k-means is different: Hennig (2007),
bootstrap Jaccard stability, implemented as `fpc::clusterboot`. If you want a stability
number, use that and call it stability.

What this package reports instead is what k-means can honestly support: the partition,
its silhouette, what distinguishes each cluster, and how the clusters relate to
variables you measured independently.

## One cohort

```bash
kmeans-py choose-k --project cohortA --k-max 8 --matrix cohortA.csv --cluster rows \
    --log2 --adjust combat --batch-col Batch --protect Group --metadata meta.csv
```

```bash
kmeans-py cluster --project cohortA --k 4 --k-rows 4 --k-cols 2 \
    --matrix cohortA.csv --cluster rows \
    --log2 --adjust combat --batch-col Batch --protect Group --metadata meta.csv \
    --annotate "Group,Sex,Age_group,Disease_activity,SLEDAI_2K" --numeric SLEDAI_2K \
    --dist minkowski --linkage ward.D2 --plot
```

## The figure

`--plot` draws the two-way split heatmap. **Both axes are always clustered.** k-means
decides the blocks, `--k-rows` on the samples and `--k-cols` on the objects;
hierarchical clustering decides the order inside each block and draws that block its
own small dendrogram. White rules separate the blocks, and the demographics you name
in `--annotate` are tiled down the side with a colour key.

The default pairing is minkowski with ward.D2, on both axes. In pvclust, minkowski
means p=2, that is euclidean, because `dist.pvclust` never forwards `p` to R's `dist`.

A block dendrogram is local to its block. k-means is a partition, so the blocks carry
no nesting and no branch lengths between them; those trees say nothing about how one
block relates to another.

Writes the partition, the assignment, the centroids, the per-cluster composition, the
cluster-variable associations, and the objects that most distinguish each cluster.

## Reading the output

**Silhouette before anything else.** It is how much closer a member sits to its own
cluster than to the next nearest. Near 1 is well separated; near 0 means the partition
is cutting through a continuum, which is what most biological data looks like. A low
silhouette does not mean the clusters are useless, but it does mean they are a
convenient summary rather than discovered groups.

**The associations say what the split tracks, not whether it is real.** The clusters
are taken as given. And if a variable drove the clustering, finding it associated
afterwards is circular: read it as "this partition is largely that effect".

**Watch for genotype.** On the SLE cohort at k=4, the smallest cluster was defined
almost entirely by five SIGLEC5 and SIGLEC14 reagents sitting a hundred-fold below
everyone else, with tight spread on both sides. That is the SIGLEC14 null
polymorphism, not disease biology. k-means will happily hand you a clean cluster that
is a genotype, a batch, or a plate.

## Federating

The partition federates exactly, and this part holds regardless of what you think about
*p*-values. Lloyd's algorithm averages points into centroids, so it needs coordinates,
and the pooled coordinates would be every project's samples side by side. The objective
k-means minimises depends only on pairwise distances (Huygens' theorem), and
`pvclust-py` recovers the pooled distance matrix exactly from sufficient statistics. So
the partition is found centrally with no coordinates existing anywhere.

```bash
pvclust-py project-stats --project cohortA --matrix cohortA.csv ...
kmeans-py  aggregate --labels cohortA_labels.txt --k 10 --stats cohortA_stats.npz ...
```

**The privacy rule is `p < n`.** The Gram matrix has rank `min(n, p)`, so a project with
fewer samples than objects exposes its row space completely. At `n = 1` the Gram is
rank one and returns the row itself (Homer et al., 2008).

## Licence

GPL-3.0-or-later, inherited from R `pvclust` via `pvclust-py`.

## References

Hennig, C. (2007). Cluster-wise assessment of cluster stability.
*Computational Statistics & Data Analysis* 52(1):258–271.

Johnson, W.E., Li, C. & Rabinovic, A. (2007). Adjusting batch effects in microarray
expression data using empirical Bayes methods. *Biostatistics* 8(1):118–127.

Rousseeuw, P.J. (1987). Silhouettes: a graphical aid to the interpretation and
validation of cluster analysis. *J. Comput. Appl. Math.* 20:53–65.

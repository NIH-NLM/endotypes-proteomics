"""
Command line for kmeans-py.

Deliberately narrow. Reading a matrix, batch correction, the shared feature
vocabulary, sufficient statistics and the pooled distance matrix all belong to
``pvclust-py`` and are called from there, so the two packages never hold two versions
of the same thing.

    kmeans-py choose-k    how separated the clusters are at each k
    kmeans-py cluster     the partition, its profile, and the demographics
    kmeans-py aggregate   the pooled partition, from distances alone
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer

from pvclust_py.cli import (_ADJREPORT, _ADJUST, _ANNOTATE, _BATCHCOL, _DIST,
                            _FEATKEY, _FEATLABEL, _FEATMAP, _LINK, _LOG2, _MATRIX,
                            _METADATA, _PROTECT, _RFU, _SAMPLES, _SEED, _SHARED,
                            _SOMAMERS, _TOPVAR, apply_adjust, load_inputs)

app = typer.Typer(add_completion=False, help=__doc__)

_CLUSTER = typer.Option(
    "rows", "--cluster",
    help="Which axis to cluster: rows (samples -- what you want when the clusters "
         "are to be described by demographics) or columns (the measured objects)")
_K = typer.Option(..., "--k", help="Number of clusters")
_KROWS = typer.Option(None, "--k-rows",
    help="Blocks to split the SAMPLES into, for the two-way heatmap. Defaults to --k")
_KCOLS = typer.Option(None, "--k-cols",
    help="Blocks to split the OBJECTS into, for the two-way heatmap. Defaults to 2")
_DISTK = typer.Option("minkowski", "--dist",
    help="Distance for the ordering WITHIN each block. minkowski means p=2, i.e. "
         "euclidean, matching what pvclust's dist.pvclust actually does")
_LINKK = typer.Option("ward.D2", "--linkage",
    help="Agglomeration for the ordering within each block")
_NINIT = typer.Option(10, "--n-init",
                      help="k-means restarts. Lowering it makes the answer depend on "
                           "the seed more than on the data")


def _annotations(metadata, rfu, samples, somamers, feature_label):
    if metadata:
        from pvclust_py.io import read_metadata
        return read_metadata(metadata)
    if samples:
        from pvclust_py.somascan import read_somascan
        return read_somascan(rfu, samples, somamers,
                             somamer_id=feature_label or "SeqId")[1]
    return None


def _prepare(matrix, rfu, samples, somamers, top_variable, shared_features, log2,
             feature_map, feature_key, feature_label, metadata, adjust, batch_col,
             protect, adjust_report):
    """Load, rename, restrict and batch-correct -- all of it pvclust-py's code."""
    X = load_inputs(matrix, rfu, samples, somamers, top_variable, "columns",
                    log2, feature_map, feature_key, feature_label, shared_features)
    ann = _annotations(metadata, rfu, samples, somamers, feature_label)
    X = apply_adjust(X, ann, adjust, batch_col, protect, adjust_report)
    return X, (ann.reindex(X.index) if ann is not None else None)


@app.command("choose-k")
def choose_k_command(
    project: str = typer.Option(..., "--project", help="Project/cohort id"),
    k_max: int = typer.Option(10, "--k-max", help="Try every k from 2 to this"),
    matrix: Optional[Path] = _MATRIX,
    rfu: Optional[Path] = _RFU,
    samples: Optional[Path] = _SAMPLES,
    somamers: Optional[Path] = _SOMAMERS,
    cluster: str = _CLUSTER,
    top_variable: Optional[int] = _TOPVAR,
    shared_features: Optional[Path] = _SHARED,
    metadata: Optional[Path] = _METADATA,
    log2: bool = _LOG2,
    feature_map: Optional[Path] = _FEATMAP,
    feature_key: Optional[str] = _FEATKEY,
    feature_label: Optional[str] = _FEATLABEL,
    adjust: str = _ADJUST,
    batch_col: Optional[str] = _BATCHCOL,
    protect: Optional[str] = _PROTECT,
    adjust_report: bool = _ADJREPORT,
    seed: int = _SEED,
    n_init: int = _NINIT,
):
    """Silhouette and inertia at each k, so the choice is made in the open.

    No k is "significant". Silhouette peaks where the clusters are best separated;
    inertia always falls with k, so it only tells you where the fall slows. Read both,
    then say which k you used and why.
    """
    from .core import choose_k

    X, _ann = _prepare(matrix, rfu, samples, somamers, top_variable, shared_features,
                       log2, feature_map, feature_key, feature_label, metadata,
                       adjust, batch_col, protect, adjust_report)
    out = choose_k(X, range(2, k_max + 1), cluster=cluster, seed=seed, n_init=n_init)
    out.to_csv(f"{project}_choose_k.csv", index=False)
    typer.echo(out.to_string(index=False))
    best = out.loc[out["silhouette"].idxmax()]
    typer.echo(f"\nbest separated: k={int(best['k'])} "
               f"(silhouette {best['silhouette']:.3f}) -> {project}_choose_k.csv")


@app.command("cluster")
def cluster_command(
    project: str = typer.Option(..., "--project", help="Project/cohort id"),
    k: int = _K,
    matrix: Optional[Path] = _MATRIX,
    rfu: Optional[Path] = _RFU,
    samples: Optional[Path] = _SAMPLES,
    somamers: Optional[Path] = _SOMAMERS,
    cluster: str = _CLUSTER,
    top_variable: Optional[int] = _TOPVAR,
    shared_features: Optional[Path] = _SHARED,
    metadata: Optional[Path] = _METADATA,
    log2: bool = _LOG2,
    feature_map: Optional[Path] = _FEATMAP,
    feature_key: Optional[str] = _FEATKEY,
    feature_label: Optional[str] = _FEATLABEL,
    adjust: str = _ADJUST,
    batch_col: Optional[str] = _BATCHCOL,
    protect: Optional[str] = _PROTECT,
    adjust_report: bool = _ADJREPORT,
    annotate: Optional[str] = _ANNOTATE,
    numeric: Optional[str] = typer.Option(
        None, "--numeric",
        help="Comma-separated metadata columns to treat as continuous, e.g. "
             "'SLEDAI_2K,Age'. Inferred from dtype when omitted"),
    top: int = typer.Option(10, "--top",
                            help="Objects to report per cluster as distinguishing it"),
    seed: int = _SEED,
    n_init: int = _NINIT,
    plot: bool = typer.Option(False, "--plot", help="Draw the two-way heatmap"),
    k_rows: Optional[int] = _KROWS,
    k_cols: Optional[int] = _KCOLS,
    heat_dist: str = _DISTK,
    heat_linkage: str = _LINKK,
    max_labels: int = typer.Option(60, "--max-labels",
        help="Draw tick labels while an axis has at most this many objects"),
):
    """Partition the data, then describe the clusters using the demographics.

    Writes the partition, the per-cluster composition, the cluster-variable
    associations, and the objects that most distinguish each cluster.

    The associations test whether the split TRACKS a variable. They say nothing about
    whether the cluster is real -- the clusters are taken as given. And if a variable
    drove the clustering, finding it associated afterwards is circular: read it as
    "this partition is largely that effect", not as confirmation.
    """
    from .core import kmeans
    from .profile import associations, composition, distinguishing

    X, ann = _prepare(matrix, rfu, samples, somamers, top_variable, shared_features,
                      log2, feature_map, feature_key, feature_label, metadata,
                      adjust, batch_col, protect, adjust_report)
    res = kmeans(X, k, cluster=cluster, seed=seed, n_init=n_init)

    res.table().to_csv(f"{project}_clusters.csv", index=False)
    res.assignment.to_csv(f"{project}_assignment.csv")
    res.centroids.to_csv(f"{project}_centroids.csv")
    typer.echo(f"{project}: k={k} over the {cluster}, "
               f"overall silhouette {res.silhouette['overall']:.3f}")
    typer.echo(res.table()[["cluster", "n_members", "silhouette"]].to_string(index=False))

    dist = distinguishing(X if cluster == "rows" else X.T, res.assignment, top=top)
    dist.to_csv(f"{project}_distinguishing.csv", index=False)

    cols = [c.strip() for c in annotate.split(",")] if annotate else None
    num = [c.strip() for c in numeric.split(",")] if numeric else None
    if ann is not None and cols:
        missing = [c for c in cols if c not in ann.columns]
        if missing:
            raise typer.BadParameter(
                f"--annotate names {missing}, which the metadata does not have; "
                f"available: {list(ann.columns)[:10]}")
        # The same subset feeds the figure. Without this the strips show every
        # metadata column, including 116-level donor ids, which is unreadable.
        ann = ann[cols]
    if ann is not None and cluster == "rows":
        use = ann
        composition(res.assignment, use, numeric=num).to_csv(
            f"{project}_composition.csv", index=False)
        assoc = associations(res.assignment, use, numeric=num)
        assoc.to_csv(f"{project}_associations.csv", index=False)
        typer.echo("\nhow the partition lines up with what you already know:")
        typer.echo(assoc.to_string(index=False))
    elif cluster == "rows":
        typer.echo("  no --metadata given, so no demographic profile was written")
    else:
        typer.echo("  clustering the columns, so the sample demographics do not apply; "
                   "pass --cluster rows to describe sample clusters")

    if plot:
        # Both axes, always. k-means gives the blocks; hierarchical clustering orders
        # within each block and draws its dendrogram, which is what the R version did.
        from .plot import split_heatmap

        kr = k_rows or k
        kc = k_cols or 2
        rows = res if cluster == "rows" else kmeans(X, kr, cluster="rows", seed=seed,
                                                    n_init=n_init, centroid_tree=False)
        cols = res if cluster == "columns" else kmeans(X, kc, cluster="columns",
                                                       seed=seed, n_init=n_init,
                                                       centroid_tree=False)
        split_heatmap(X, f"{project}_heatmap_kmeans",
                      row_assign=rows.assignment, col_assign=cols.assignment,
                      row_annotations=ann, method_dist=heat_dist,
                      method_hclust=heat_linkage, max_labels=max_labels,
                      title=f"{project} — k-means blocks "
                            f"({rows.k} x {cols.k}), {heat_dist} + {heat_linkage}")
        typer.echo(f"  wrote {project}_heatmap_kmeans.(png|svg|html) "
                   f"— {rows.k} row blocks x {cols.k} column blocks")


@app.command("aggregate")
def aggregate_command(
    labels: Path = typer.Option(..., "--labels",
        help="Object labels, one per line (from `pvclust-py project-stats`)"),
    k: int = _K,
    stats: List[Path] = typer.Option(..., "--stats",
        help="Per-project stats npz from `pvclust-py project-stats` (repeat)"),
    dist: str = _DIST,
    linkage: str = _LINK,
    seed: int = _SEED,
):
    """Pool the projects into one partition, from distances alone.

    k-medoids, not Lloyd's. Lloyd's averages points into centroids, so it needs
    coordinates -- and the pooled coordinates would be every project's samples side by
    side, which is the one thing that must not leave a project. The objective k-means
    minimises depends only on pairwise distances (Huygens' theorem), and pvclust-py
    recovers the pooled distance matrix exactly from sufficient statistics. So the
    partition is found centrally without the coordinates existing anywhere.
    """
    import pandas as pd

    from pvclust_py.aggregate import federated_tree
    from pvclust_py.io import read_stats

    from .aggregate import federated_kmeans

    names = [l for l in Path(labels).read_text().splitlines() if l]
    D, _Z, _tree = federated_tree([read_stats(s) for s in stats], names,
                                  method_dist=dist, method_hclust=linkage)
    catalogue = federated_kmeans(D, names, k=k, seed=seed)
    pd.DataFrame(D, index=names, columns=names).to_csv("federated_distance.csv")
    pd.DataFrame([{**e, "members": ";".join(e["members"])} for e in catalogue]
                 ).to_csv("federated_kmeans_partition.csv", index=False)
    typer.echo(f"pooled {len(stats)} projects, k={k} -> federated_distance.csv, "
               f"federated_kmeans_partition.csv ({len(catalogue)} clusters)")


@app.command("project-profiles")
def project_profiles_command(
    project: str = typer.Option(..., "--project", help="Project/cohort id"),
    k: int = _K,
    matrix: Optional[Path] = _MATRIX,
    rfu: Optional[Path] = _RFU,
    samples: Optional[Path] = _SAMPLES,
    somamers: Optional[Path] = _SOMAMERS,
    shared_features: Optional[Path] = _SHARED,
    metadata: Optional[Path] = _METADATA,
    log2: bool = _LOG2,
    feature_map: Optional[Path] = _FEATMAP,
    feature_key: Optional[str] = _FEATKEY,
    feature_label: Optional[str] = _FEATLABEL,
    adjust: str = _ADJUST,
    batch_col: Optional[str] = _BATCHCOL,
    protect: Optional[str] = _PROTECT,
    adjust_report: bool = _ADJREPORT,
    annotate: Optional[str] = _ANNOTATE,
    seed: int = _SEED,
    n_init: int = _NINIT,
    min_cluster: int = typer.Option(10, "--min-cluster",
        help="Refuse to release if any cluster is smaller than this"),
):
    """Partition THIS project's patients with k-means, and write the releasable summary.

    The k-means counterpart of `pvclust-py project-profiles`, which uses pvpick and
    needs no k. Here you choose k, so choose it in the open: run `choose-k` first and
    say what you picked and why.

    Clustering patients cannot be federated the way clustering features is: projects
    hold disjoint patients, so there is no shared object vocabulary and a cluster is a
    set of subject ids. What crosses the boundary is a per-cluster MEAN PROFILE and a
    member count -- an average over patients, carrying no subject id.

    Two release rules apply. Everything released from one dataset is ONE release: ship
    the k you chose, not every k you tried, because two partitions of the same patients
    differing by one member reveal that member by subtraction. And a partition with any
    cluster below --min-cluster is suppressed WHOLE, because dropping only the small
    cluster announces it existed.
    """
    import pandas as pd

    from pvclust_py.profiles import SuppressedError, project_profiles

    from .core import kmeans

    X = load_inputs(matrix, rfu, samples, somamers, None, "rows",
                    log2, feature_map, feature_key, feature_label, shared_features)
    ann = _annotations(metadata, rfu, samples, somamers, feature_label)
    X = apply_adjust(X, ann, adjust, batch_col, protect, adjust_report)
    if ann is not None:
        ann = ann.reindex(X.index)
        if annotate:
            ann = ann[[c.strip() for c in annotate.split(",")]]

    res = kmeans(X, k, cluster="rows", seed=seed, n_init=n_init, centroid_tree=False)
    try:
        out = project_profiles(X, res.assignment, project=project,
                               min_cluster=min_cluster, annotations=ann)
    except SuppressedError as e:
        raise typer.BadParameter(str(e))

    out["profiles"].to_csv(f"{project}_profiles.csv")
    out["counts"].to_csv(f"{project}_profile_counts.csv", index=False)
    res.assignment.rename("cluster").to_csv(f"{project}_patient_assignment.csv")
    if out["demographics"] is not None:
        out["demographics"].to_csv(f"{project}_profile_demographics.csv", index=False)

    typer.echo(f"{project}: k={k} over {X.shape[0]} patients, "
               f"sizes {list(out['counts']['n_patients'])}, "
               f"silhouette {res.silhouette['overall']:.3f}")
    typer.echo(f"  wrote {project}_profiles.csv -- no subject leaves")


def main():
    app()


if __name__ == "__main__":
    main()

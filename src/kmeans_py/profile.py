"""
Describing the clusters with the demographics you already hold.

This is what replaces a *p*-value on the cluster itself. k-means will always return k
clusters, and nothing in the algorithm says whether they mean anything. What can be
said honestly is whether they line up with something known -- disease group, sex, age
band, disease activity -- and how strongly.

Two views, because they answer different questions:

  COMPOSITION  what is in each cluster. A plain cross-tabulation, as counts and as
               row percentages. This is what you read to describe a cluster.

  ASSOCIATION  whether the split tracks a variable at all, across every cluster at
               once, with an effect size and a p-value. Note what is being tested:
               the *p*-value is about the cluster-variable association, NOT about
               whether the cluster is real. The clusters are taken as given.

A caution that matters when the clusters were built from the same data. If a variable
drove the clustering, finding it associated afterwards is circular. Read a strong
association as "this partition is largely the sex effect", not as confirmation.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd


def composition(assignment: pd.Series, annotations: pd.DataFrame, *,
                numeric: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """What each cluster contains, one row per cluster and variable level.

    Categorical variables give counts and the within-cluster percentage. Numeric ones
    give the median and the interquartile range, since clinical scores are usually
    skewed and a mean would misdescribe them.
    """
    ann = annotations.reindex(assignment.index)
    numeric = set(numeric or [c for c in ann.columns
                              if pd.api.types.is_numeric_dtype(ann[c])])
    rows = []
    for cid, idx in assignment.groupby(assignment).groups.items():
        sub = ann.loc[idx]
        for col in ann.columns:
            if col in numeric:
                v = pd.to_numeric(sub[col], errors="coerce").dropna()
                if v.empty:
                    continue
                rows.append({"cluster": str(cid), "variable": col, "level": "",
                             "n": int(v.size), "pct": np.nan,
                             "median": float(v.median()),
                             "iqr_low": float(v.quantile(0.25)),
                             "iqr_high": float(v.quantile(0.75))})
            else:
                counts = sub[col].astype(str).value_counts()
                for level, n in counts.items():
                    rows.append({"cluster": str(cid), "variable": col, "level": level,
                                 "n": int(n), "pct": 100.0 * n / len(sub),
                                 "median": np.nan, "iqr_low": np.nan,
                                 "iqr_high": np.nan})
    return pd.DataFrame(rows)


def associations(assignment: pd.Series, annotations: pd.DataFrame, *,
                 numeric: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Does the partition track each variable? Effect size and *p*-value per variable.

    Delegates to :func:`pvclust_py.diagnostics.association`, so the two packages score
    a cluster-variable association identically and the numbers can be compared. Cramer's
    V is bias-corrected (Bergsma 2013): the uncorrected form inflates badly when a
    variable has many levels relative to the sample size, and scored a 96-level plate
    position at 0.56 on data whose chi-square *p*-value was 0.46.
    """
    from pvclust_py.diagnostics import association

    return association(assignment, annotations.reindex(assignment.index),
                       numeric=numeric)


def distinguishing(X, assignment: pd.Series, *, top: int = 10) -> pd.DataFrame:
    """Which measured objects separate each cluster from the rest.

    A standardised mean difference per cluster and object: the cluster's mean minus
    everyone else's, over the pooled standard deviation. Positive means raised in that
    cluster. This is descriptive and deliberately has no *p*-value attached -- every
    object was used to build the partition, so testing it against that same partition
    would be circular.
    """
    if not isinstance(X, pd.DataFrame):
        X = pd.DataFrame(X)
    A = X.reindex(assignment.index)
    rows = []
    for cid in sorted(assignment.unique()):
        inside = A[assignment == cid]
        outside = A[assignment != cid]
        if inside.empty or outside.empty:
            continue
        pooled = np.sqrt((inside.var(ddof=1).fillna(0) +
                          outside.var(ddof=1).fillna(0)) / 2).replace(0, np.nan)
        d = ((inside.mean() - outside.mean()) / pooled).dropna()
        for obj, val in d.reindex(d.abs().sort_values(ascending=False).index[:top]).items():
            rows.append({"cluster": str(cid), "object": obj,
                         "std_mean_difference": float(val),
                         "direction": "higher" if val > 0 else "lower"})
    return pd.DataFrame(rows)

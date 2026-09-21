# endotypes-proteomics

**Unsupervised protein modules, then patient endotypes, in plasma proteomics — one step per
notebook, one step per container.**

> The directory on disk is still `kmeans-py`, from an earlier phase when this was a k-means
> experiment in Python. Neither half is accurate any more. See [Naming](#naming).

## What this does

Given a protein abundance matrix and patient metadata, it asks whether there are **endotypes** —
sets of patients sharing a coordinated shift across a set of proteins. That is a *block*, a
patient × protein rectangle, not two dendrograms glued to the sides of a heatmap.

The order is the argument:

1. **Define protein modules first, unsupervised**, told nothing about the patients.
2. **Then** summarise each module as one number per patient and look for patient structure.

Selecting proteins by how well they separate known groups collapses them onto one axis — 60 proteins
selected that way correlated at median 0.73 and one component explained 72% of their variance, so
there were no modules left to find. Defining modules first avoids that.

## The data

SomaScan v4.1 plasma proteomics in systemic lupus erythematosus, Imperial College London
(Leung, Pickering, Botto, Peters).

> **https://doi.org/10.5281/zenodo.20342569**

Download the record and unzip it so the files sit at `data/SLE_doi.10.5281_zenodo_20342569/`:

```
data/SLE_doi.10.5281_zenodo_20342569/
├── abundance.csv          369 samples x 7,288 SOMAmers, RFU
├── sample-metadata.csv    Group (283 SLE / 86 HV), Included_in_study, Batch,
│                          Disease_activity, SLEDAI_2K, C3/C4, autoantibody status
└── feature_metadata.txt   SeqId -> Target / UniProt / GeneSymbol
```

`data/` is gitignored, and so is everything in `ipynb/` derived from it. **No study data is
committed.** The one exception is `ipynb/cohort_assignment.csv`, which freezes cohort membership and
carries sample identifiers but no abundance values.

## Install and run

```bash
micromamba env create -f endotypes-proteomics.yml
micromamba activate endotypes-proteomics
Rscript -e 'IRkernel::installspec()'    # R kernel
python -m bash_kernel.install           # Bash kernel
```

By hand:

```bash
jupyter lab
```

By machine:

```bash
./run_all.sh              # steps 00-07
./run_all.sh --optional   # also 08 and 09, the projections
```

All analysis is R. JupyterLab, `nbconvert` and `bash_kernel` are Python packages, so conda pulls a
Python interpreter in as their dependency — it is the notebook host, not a language this project
writes code in.

## The steps

Each notebook is one stage, reads the artifact the previous stage wrote, and writes its own. No
stage refits what an earlier stage fitted, which is what makes them separable, cacheable and
containerizable one per stage.

| notebook | does |
|---|---|
| `00_prepare_data` | Zenodo → filter → log2 → ComBat → annotate → de-duplicate → split |
| `01_soft_threshold` | the soft-thresholding power |
| `02_modules` | WGCNA fit, unsupervised |
| `03_eigenproteins` | one number per patient per module |
| `04_module_traits` | module ↔ clinical trait, BH-corrected |
| `05_endotypes` | cluster patients in module space |
| `06_heatmap` | the patient × protein figure |
| `07_federation` | what would cross an institutional boundary |
| `08_project_healthy` | the 86 healthy volunteers, scored on SLE-defined modules |
| `09_project_timepoints` | *optional* — later visits of repeat donors |

### Correct once, before anything

Step 00 does **every** fix to the data in one pass, across all 356 included samples, before any
biological selection and before any split. Correction is a property of the assay, not of the
question. In order: log2 → ComBat → gene symbols → age as an ordered band → one sample per donor →
split.

**Healthy volunteers stay in for the correction and are then held out of every fit.** They are
*projected* in step 08, not pooled — see below.

**One sample per donor.** 270 SLE samples come from 260 donors; six were sampled repeatedly and
`SampleId` encodes the visit. The first visit is kept, and it is also the most active disease
(D256: SLEDAI-2K 9 → 6 → 4 → 2; D173: 12 → 5), so "first by date" and "highest activity" pick the
same sample. Later visits go to their own matrix for step 09.

**Cohort membership is frozen** to `cohort_assignment.csv` and read back on every run, so no
downstream number changes because of a reshuffle. Cohorts are 87 / 87 / 86, stratified within batch.

## What was found

### The published interferon module is recovered

`ivory`, 14 proteins, 10 of them canonical interferon-stimulated: STAT1 (×2), DDX58, IFIT3, ISG15
(×2), GBP1, MX1, IFIH1, CXCL10. Hub by module membership is **ISG15 `seq.14151.4`** at kME 0.94 —
ISG15 is the hub protein the source paper names for its own interferon module, reached on n = 207 by
a different route.

It is also the most clinically connected module in the panel — nine associations surviving BH
correction across 780 tests:

| trait | r | q |
|---|---|---|
| **anti-Sm** | +0.53 | 7.5e-05 |
| C3 | −0.45 | 1.5e-03 |
| age band | −0.45 | 1.6e-03 |
| **anti-Ro60** | +0.40 | 8.7e-03 |
| disease duration | −0.39 | 8.8e-03 |
| anti-dsDNA | +0.39 | 1.2e-02 |
| **anti-RNP-A** | +0.37 | 1.6e-02 |
| **anti-RNP68** | +0.36 | 2.6e-02 |
| lymphocyte count | −0.35 | 3.4e-02 |

The paper reports that antibodies to **RNA-binding proteins — anti-Sm, anti-Ro60, anti-RNP68,
anti-RNP-A** — go with raised interferon-stimulated proteins. All four are here and all four survive
correction; anti-Sm is the strongest association in the whole table. **Anti-La, which is not on that
list, is flat** (r = +0.16, q = 0.65) — the one antibody that should not track is the one that does
not.

### Healthy volunteers, projected rather than pooled

WGCNA defines modules by correlation **across samples**, so pooling healthy volunteers in lets the
between-group mean difference manufacture correlation: any two proteins elevated in disease
correlate whether or not they are co-regulated, and you get a large module that is really the group
contrast. Step 08 instead scores the healthy volunteers on axes they had no part in defining.

The interferon module separates them at **−1.02 SLE standard deviations, q = 1.9e-12** — the
strongest specific separation in the panel, and it does it on 14 proteins. The overlap is the
interesting part: **4 of 86 healthy volunteers sit above the SLE median, and 10 of 87 SLE patients
sit below the healthy median.** The continuum, counted.

The projection is guarded: pushing the SLE matrix back through the same path must reproduce
`moduleEigengenes`, and it does at **r = 1.0000**. Without that gate the code could be refitting
rather than projecting, and every number after it would be meaningless.

**This separation is a lower bound.** ComBat ran unsupervised (`mod = NULL`), and healthy volunteers
are 17% of batch A but 44% of batch B, so part of the true batch shift is group composition and the
correction absorbs some of the healthy-vs-SLE difference with it. Protecting `Group` instead would
make any later unsupervised separation partly an artifact of having protected it. Under-detecting is
the right direction to err.

### The patient partition is weak, and interferon is not in it

Best silhouette 0.20 at k = 2 — below the ~0.25 convention for substantial structure. The split is
driven by `blue` (1,213 proteins) and `black` (150), i.e. one dominant axis. **The interferon module
separates the two endotypes by about 0.1 SD, essentially not at all**, despite carrying the
strongest clinical associations. The interferon axis is orthogonal to the partition, so those
patients are a group this clustering does not isolate. Clustering on the interferon module alone is
the obvious next move.

Guarded too: an unconstrained search preferred k = 4 with a cluster of 2. **A cluster of one has
silhouette 1 by construction**, so step 05 requires the smallest cluster to hold at least ten
patients.

### Federation is described, not performed

Step 07 proves that pooling four sufficient statistics — `N`, `S`, `Q`, `G = XᵀX` — reproduces the
pooled correlation matrix to **2.5e-12**, so WGCNA federates exactly in principle. It also states
the blocker: `rank(G) = min(n, p)`, and with 7,288 proteins against ~87 patients per site, p ≫ n and
the Gram matrix can be inverted back toward the rows. The panel must be reduced below n before any
Gram leaves a site. Everything here used pooled raw data and is a **simulation of what federation
would produce**, not a federated run.

## Open work, in priority order

1. **`modulePreservation` across A/B/C.** The module's *boundary* moves with the sample: earlier
   draws of ~90 patients from the same cohort put these proteins in modules of 12, 23, 39 and once
   775 members, and the module count at identical parameters has ranged 25–52. The proteins
   co-cluster reliably; where the boundary falls does not. Until this runs, the table above
   describes cohort A.
2. **Support layer**: bootstrap Jaccard per patient cluster (`fpc::clusterboot`); permutation
   p-value per patient × protein block. Silhouette measures *separation*; neither it nor a q-value
   measures whether a cluster would reappear.
3. **VarSelLCM** replacing the hardcoded 15-trait list in step 04, with the snRNP group held out —
   anti-Sm correlates 0.50 with anti-RNP-A and 0.44 with anti-RNP68, so holding out anti-Sm alone
   leaks.
4. **Tiered WGCNA** on preserved modules above 200 proteins, one pass, no recursion.
5. **Cluster on the interferon module alone**, since it is orthogonal to the current partition.
6. **Finish the plasma-protein annotation** (5,401 of 6,401 reagents; the missing ones include ALB,
   CRP, C3, SAA1, HP and ISG15). Modules hubbed on implausible proteins — ACRV1, a sperm acrosomal
   protein, topped one — cannot be filtered without it.

## Prior art: WGCNA on SomaScan

- **Leung *et al.*, plasma proteomic architecture of SLE** — the source of this dataset. SomaScan
  v4.1, 7,288 analytes, 21 protein modules on n = 207 SLE, interferon module hub ISG15.
  [JCI Insight](https://insight.jci.org/articles/view/206938) ·
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13461151/)
- **CSF proteome in frontotemporal lobar degeneration** — aptamer proteomics, 31 modules.
  [Nature Aging](https://www.nature.com/articles/s43587-025-00878-2)
- **Protein-specific co-expression in Alzheimer's disease.**
  [Cell Systems](https://www.cell.com/fulltext/S2405-4712(16)30370-2)
- **Reproducible molecular subtypes in SLE.**
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13332139/)

All work at n in the low hundreds or more. None attempts module discovery at n = 87.

## Naming

**`endotypes-proteomics`** carries no method and no language, deliberately. Two methods are in play —
WGCNA for the proteins, VarSelLCM for the clinical variables — so a method name would have to pick
between them. A method name is also a hostage: WGCNA is **Langfelder & Horvath's** R package, and a
repository called plainly `wgcna` under an institutional organisation would read as though it *is*
that package. And a language suffix is a promise the sibling Nextflow repositories would break —
exactly how `kmeans-py` became wrong.

The k-means implementation (`src/kmeans_py/`, `tests/`, `scripts/`, `pyproject.toml`, `ipynb/`) has
been deleted. It remains in git history, which is the point of renaming rather than starting a new
repository: the record of what was tried and rejected is part of the argument for what replaced it.

Sibling repositories: `pvclust-py` (hierarchical clustering with AU *p*-values) and `oadr-cpep`
(federated supervised regression).

## Why `pvclust` is not used here

`pvclust` attaches an approximately-unbiased *p*-value to each cluster of a dendrogram. It is the
right tool for *"is this cluster real"* and the wrong tool for this question:

1. **It answers about one axis at a time.** Two trees run separately do not produce a block.
2. **On this orientation its p-values are anti-conservative.** Clustering patients resamples
   proteins, and proteins are co-expressed rather than independent draws. On the patient profiles
   every edge returned AU = 1.000 with `df = 0` — the curve was never fitted and the number printed
   is not a p-value.
3. **It never chooses the proteins.** Feeding it the most-variable proteins made the answer a
   consequence of that choice, and the most-variable proteins in plasma are the abundant ones.

It remains the right tool for assessing the **protein** tree, and lives in
[`pvclust-py`](https://github.com/NIH-NLM/pvclust-py).

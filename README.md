# endotypes-proteomics

**Unsupervised protein modules, then patient endotypes, in plasma proteomics — one step per
notebook, one step per container.**

Given a protein abundance matrix and patient metadata, we are looking for  **endotypes** —
patients sharing a similar pattern of protein expression. 

## Conceptual approach
Proteins are co-expressed additionally with a measurement technology such as the SomaScan array
There are multiple probes to the same protein.  

1. **Define protein modules unsupervised** -- Protein module clustering without patient information.
2. **Describe protein modulesn** Look for patterns in patient structure within these protein modules.
   
## Example

SomaScan v4.1 plasma proteomics in systemic lupus erythematosus, Imperial College London
(Leung, Pickering, Botto, Peters).

**https://doi.org/10.5281/zenodo.20342569**

## Prepare the data

### Running on ADAPTS

1. Login to lifebit.ai
2. Use your PIV card to authenticate
3. Spin up JupyterLab Notebook (8 vCPUs, 16GB RAM)

When ready, open a terminal window

4. Clone the repository and change directory 

```bash
git clone https://github.com/NIH-NLM/endotypes-proteomics.git
cd endotypes-proteomics
```

5. Create a directory under the data directory

```bash
mkdir -p data/SLE_doi.10.5281_zenodo_20342569
cd data/SLE_doi.10.5281_zenodo_20342569
curl https://zenodo.org/api/records/20342569/files-archive > SLE_doi.10.5281_zenodo_2034256.zip
unzip SLE_doi.10.5281_zenodo_2034256.zip
rm SLE_doi.10.5281_zenodo_2034256.zip
cd ../..
pwd
```

6. Confirm you have the data.
```bash
ls -l data
```
The files should now show these three files

```bash
data/SLE_doi.10.5281_zenodo_20342569/
├── abundance.csv          369 samples x 7,288 SOMAmers, RFU
├── sample-metadata.csv    Group (283 SLE / 86 HV), Included_in_study, Batch,
│                          Disease_activity, SLEDAI_2K, C3/C4, autoantibody status
└── feature_metadata.txt   SeqId -> Target / UniProt / GeneSymbol
```

### Repository layout

Two kinds of file, and the split is the point: what a run **produces** is disposable, what a run
**depends on** is committed.

```
cohorts/     COMMITTED  the frozen splits, cohort_assignment.csv, clinical-traits.csv
proteins/    COMMITTED  curated prior knowledge -- interferon-response-genes.json
ipynb/       COMMITTED  the ten notebooks, and nothing else
src/         COMMITTED  R that is not a notebook, including paths.R
data/        ignored    the Zenodo download
  run_artifacts/        every .rds, derived .csv and figure a run makes -- delete freely
```

`cohorts/` is committed because redrawing the split would move every number downstream of it;
`data/run_artifacts/` is not, because `./run_all.sh` rebuilds all of it from `cohorts/`.

**`src/paths.R` is the only place these locations are written down.** Every notebook opens with
`source("../src/paths.R")` and then says `coh("R_cohort-A_meta.csv")` or `art("wgcna_A.rds")`
rather than a relative path. Moving a directory is a one-line change there.

### Environment Setup

```bash
mamba env create -f endotypes-proteomics.yml     # or micromamba / conda
mamba activate endotypes-proteomics
Rscript -e 'IRkernel::installspec()'             # R kernel
python -m bash_kernel.install                    # Bash kernel
jupyter lab
```

**Two R packages are not in the environment file, deliberately.** `WGCNA` has no
`osx-arm64` conda build — declaring it makes the environment unsolvable on Apple Silicon — and
`VarSelLCM` has no conda package on any platform. Both are installed from CRAN by `ensure_pkg()`
in `src/paths.R`, which every notebook that needs them calls before `library()`. Nothing extra to
run: open any notebook under the R kernel and it resolves its own dependencies, installing only
what is missing and touching no network on a second run.

**One dependency is not in that file.** `VarSelLCM` has no conda package, so step 04 installs it
from CRAN in its own first cell, guarded by `requireNamespace()` so a second run is a no-op. It is
the only thing this project installs from inside a notebook, and the environment file says so where
the dependency would otherwise have gone.

### Notebooks

For ease of analysis and understanding, a notebook was created for each of the separate execution steps
Each notebook reads the artifact of the previous step and then writes out its own output.

**The sequence runs from a clean clone.** After the download above, `./run_all.sh` executes steps
00–07 and 10–16 in order with no manual intervention; `--optional` adds 08 and 09. Verified by
deleting `data/run_artifacts/` entirely and re-running: every step completes, and **`cohorts/`
regenerates bit-identically** — step 00 reads the committed `cohort_assignment.csv` and reseeds
from it, so the split never moves.

Only `data/` is required from outside the repository. Everything else is either committed
(`cohorts/`, `proteins/`, `ipynb/`, `src/`) or regenerated into `data/run_artifacts/`.

| notebook | does |
|---|---|
| `00_prepare_data` | Zenodo → filter → log2 → ComBat → annotate → de-duplicate → split |
| `01_soft_threshold` | the soft-thresholding power |
| `02_modules` | WGCNA fit, unsupervised |
| `03_eigenproteins` | one number per patient per module |
| `04_module_traits` | module ↔ clinical trait, BH-corrected; VarSelLCM on the trait list |
| `05_endotypes` | cluster patients in module space |
| `06_heatmap` | the patient × protein figure |
| `07_federation` | what would cross an institutional boundary |
| `10_federated_modules` | performs it: one pooled definition, reapplied to A, B, C |
| `11_modules_per_cohort` | the same fit for B and C — *slow, two full WGCNA runs* |
| `12_panels_per_cohort` | per-cohort panels under three trait conditions |
| `13_cluster_both_axes` | ward.D2 / Minkowski and k-means, both axes, modules not imposed |
| `14_heatmaps` | two-tier figures: overview, then labelled zooms, shared trait legend |
| `15_kmeans_arm` | the k-means comparator, drawn against the hierarchical arm |
| `16_project_and_federate` | cohort B by projection, and the federation arithmetic |
| `17_cohort_diagnostics` | why the three cohorts differ |
| `18_federate_per_module` | federation performed, one module at a time |
| `19_federation_benefit` | does it help? three levels, measured |
| `20_module_preservation` | do the modules exist in the other cohorts? |
| `08_project_healthy` | the 86 healthy volunteers, scored on SLE-defined modules |
| `09_project_timepoints` | *optional* — later visits of repeat donors |

### Initial Data Transformation and Batch Correction

Step 00 prepares the data in one pass.

There are 6 subjects that have multiple samples.  Only the first sample is used.
Later visits are explored in step 09

### Cohorts

Three cohorts A,B and C were created by randomly sampling, so we could illustrate batch correction.
The assignment is frozen in `cohorts/cohort_assignment.csv` and read back on every run; delete that
file to redraw, deliberately.

### The two committed metadata files

**`cohorts/clinical-traits.csv`** — the 15 traits step 04 tests, one row each: `name`,
`source_column`, `type` (how to make a number of it), `group`, `description`. The `group` column is
the load-bearing one. Antibodies to the same antigen system carry the same information — in cohort A
anti-Sm correlates **+0.52** with anti-RNP-A and anti-Ro52 correlates **+0.73** with anti-Ro60 — so
a hold-out has to remove a whole group at a time. Removing anti-Sm alone leaks through anti-RNP-A.

**`proteins/interferon-response-genes.json`** — 37 curated interferon-stimulated genes, 19 of them
on this SomaScan menu carrying 23 probes, each with gene symbol, UniProt, protein name and its
SOMAmer SeqIds. It replaces `grep("14148")`, which used to appear in four notebooks: that lookup
returned whichever module one probe fell into and called it interferon, and **could not fail**.
`locate_ifn_module()` takes the module holding a plurality of the curated set and stops if no module
holds at least three of them. On cohort A it recovers `ivory` with 10 of the set against 2 for the
runner-up — a result, since the set was written from the literature rather than from this fit.

Note `feature_metadata.txt` ships `GeneSymbol` **pre-disambiguated**, so multi-probe genes arrive as
`ISG15_seq.14148.2`. Matching on a bare `"ISG15"` silently finds nothing; the JSON carries the
column name for this reason.

### Outcomes

The original authors also used WGCNA. We were able to recapitulate the findings in the original paper 
including the interferon module is recovered

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

### There are no patient endotypes in cohort A

**Zero partitions pass, at any k, in either candidate module space.** Step 05 selects k by
**prediction strength** (Tibshirani & Walther): split the patients in half, cluster each half, and
test whether one half's centroids predict co-membership in the other. The published threshold is
0.8. Cohort A scores:

| space | k = 2 | k = 3 | k = 4 | k = 5 |
|---|---|---|---|---|
| all 12 associated modules | 0.65 | 0.55 | 0.46 | 0.32 |
| 8 specific (≤ 50 proteins) | 0.69 | 0.47 | 0.40 | 0.29 |

Nothing reaches 0.8, and the score **falls monotonically** with k — data with real structure at
some k shows a peak there. Bootstrap Jaccard (Hennig) agrees: the clusters that reappear are always
the large residual group, while the small ones sit at 0.49–0.65. The figure in step 06 draws
k = 2…5 so this is visible — the big block never subdivides, it peels off 7 patients, then 2, then
2 more.

**Silhouette was removed, deliberately.** It grades k-means with k-means' own assumptions, it is not
comparable across feature spaces, and on this data it preferred a partition whose small cluster
dissolves at Jaccard 0.49. The reasoning is written out in step 05.

**The interferon module is not in the split, and that is the finding.** It separates the descriptive
groups by about 0.1 SD despite carrying the strongest clinical associations and cleanly separating
healthy volunteers in step 08. The module with the evidence behind it is the one this clustering
cannot see — which argues for clustering on the interferon module alone, not for a better k.

### Federation, performed

Step 10 runs it. Each site ships 1-D summaries only; the sites agree a panel of **p = 80 < n = 86**;
each ships `N, S, Q, G` on that panel; the pooled correlation is exact to **4.4e-12**; one module
definition is built from the correlation matrix alone and reapplied to A, B and C with federated
loadings and pooled scaling. The interferon module bands the same way at all three sites — mean
within-module r of 0.51, 0.37, 0.45 — including the two that did not define it.

**The control is the point.** On a panel chosen with no prior knowledge — the top 80 proteins by
pooled variance — only 3 of 23 ISG probes survive the reduction, and none lands in a named module.
The signature is **not recoverable**. On the ISG-anchored panel, same budget and same release rule,
it is.

> Discovery needs the full panel and therefore one site. Confirmation federates exactly, on a panel
> small enough to release, **provided the thing being confirmed was named in advance.**

That is why `proteins/interferon-response-genes.json` is committed prior knowledge rather than a
by-product of the fit: it is the object that makes the federated arm possible.

### The algebra behind it

Step 07 proves that pooling four sufficient statistics — `N`, `S`, `Q`, `G = XᵀX` — reproduces the
pooled correlation matrix to **2.5e-12**, so WGCNA federates exactly in principle. It also states
the blocker: `rank(G) = min(n, p)`, and with 7,288 proteins against ~87 patients per site, p ≫ n and
the Gram matrix can be inverted back toward the rows. The panel must be reduced below n before any
Gram leaves a site. Everything here used pooled raw data and is a **simulation of what federation
would produce**, not a federated run.

## Per-cohort discovery (steps 11–15)

Each cohort analysed separately, nothing from one selecting anything in another. Panels come from
WGCNA modules that survive a BH test against patient attributes, minus one stated size criterion:

> A module is excluded if it contains **≥1% of the assayed probes** (≥73 of 7,288), because a
> module much larger than that approximates a leading principal component of the whole panel and
> dominates the Euclidean geometry, preventing small modules from resolving.

**The threshold is chosen by a scan, not asserted.** Step 12 reports, for each candidate cut, how
well *free* hierarchical clustering on the protein axis recovers the WGCNA modules it was never
given — best Jaccard of any protein cluster against `ivory` (interferon) and `bisque4` (renal):

| panel | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 |
|---|---|---|---|---|---|---|
| <10% (765 probes) | 0.08 | 0.11 | 0.11 | 0.14 | 0.14 | **0.14** |
| <5% (392) | 0.09 | 0.15 | 0.39 | 0.39 | 0.39 | 1.00 |
| <2% (242) | 0.15 | 0.39 | 0.39 | **1.00** | 1.00 | 1.00 |
| **<1% (164)** | 0.39 | 0.39 | **1.00** | 1.00 | 1.00 | 1.00 |

At a 10% cut `ivory` is **never** recovered at any k. At 1% it is exact from k = 5 and `bisque4`
from k = 7 — which is why steps 13–15 run **k = 3…8** on both axes rather than stopping at 5.
Caveat stated plainly: choosing the panel on which free clustering reproduces the modules is a
**consistency** criterion and mildly self-confirming; the defence is the mechanism, not the
prettier answer.

| cohort | modules | trait-associated | panel at the 1% cut |
|---|---|---|---|
| A | 52 | 12 | **164 probes (144 proteins), 8 modules** |
| B | 23 | 1 (`brown`, renal, 6.7%) | **empty** — its one module exceeds the cut |
| C | 46 | 4 | 119 probes (110 proteins), 3 modules |

Step 11 refits B and C with parameters identical to step 02's cohort-A fit, so any difference
between cohorts is a difference in the data rather than in the parameters. The three fits are not
alike: B resolves into **less than half** as many modules as A and leaves 64% more protein
unassigned. Module discovery at n ≈ 87 is not stable across a random split of the same donors, and
that is the governing caveat on everything in steps 12–16.

**Cohort B has one module, and it is not a multiple-testing artifact** — cutting its tests from 345
to 138 via VarSelLCM leaves it at one. B's only clinical signal is renal.

**The modules largely survive being ignored.** Steps 13–15 cluster both axes freely, module labels
never imposed. Cohort A reaches **ARI 0.985 at k = 8** and reproduces **6 of its 8 modules exactly**
(Jaccard 1.00). `fig14_final_heatmap_A.png` is that result: every column block carries the name of
the module a free clustering found there.

**And the modules are preserved across cohorts (step 20).** `WGCNA::modulePreservation`, A as
reference, 200 permutations, read against the `gold` random-module null at Z ≈ 9:

| | in B | in C |
|---|---|---|
| preserved above the null | **7 of 8** | **7 of 8** |
| `ivory` (interferon) | Z = 10.0 | Z = 14.7, medianRank 1 |
| `brown4` | 5.9 — fails | 6.3 — fails |

**This explains why the cohorts looked so different.** A's modules *are* present in B's correlation
structure; B simply did not cut them out. The 52-vs-23 module count is a property of
`blockwiseModules` and its parameters, not of the biology.

Preservation is about proteins, significance is about patients, and step 20 never calls a coherent
protein cluster "significant".

**Projection carries structure into B, and so does federation — one module at a time.** A→B
projection recovers 5 associations, 4 of them on `bisque4`, a 10-protein renal module. The gate
holds at r = 1.0000.

Federation fails only if you try to release a whole panel as one Gram, which nothing requires.
Counted in **proteins rather than probes** and released **one module per run**, 9 of A's 11
trait-associated modules and 3 of C's 4 clear `p < n = 86` — `bisque4` is 7 proteins, `ivory` 12,
`mediumpurple3` 13. Step 17 performs it, exact to 4e-12.

**What federation buys — measured, in step 19.** Every eligible module scored in every cohort at
three levels: the cohort alone, plus the origin's **protein list**, plus **federated loadings**.

| | A | B | C | total |
|---|---|---|---|---|
| alone | 29 | **0** | 8 | 37 |
| + shared protein list | 44 | 11 | 21 | **76** |
| + federated loadings | 44 | 16 | 20 | 80 |

**The shared definition is worth everything; the federated loadings are worth nothing.** Cohort B
goes from 0 to 11 simply by being told which proteins form the module — its own fit contains no
such module. The further 76 → 80 is not a gain: the two scores correlate at **0.968–1.000**,
`|Δ|r|| > 0.05` in only **2 of 540** cells, and the traits that cross FDR 5% do so with unchanged
or *smaller* effect sizes because nine of them sit piled at q ≈ 0.046–0.052 and BH reorders them
together.

That matters for design, because the two levels have very different privacy costs. **A protein
list is a list of names.** A Gram matrix must satisfy `p < n` and be released one module at a time.
Here the privacy-cheapest option delivers essentially all of the benefit — at n ≈ 87, on modules of
7–78 proteins, where PC1 is already well determined.

**One module per release, deliberately.** Nine separate Grams totalling 222 proteins is not
obviously the same disclosure as one 222-protein Gram. Step 17 is built so that question never has
to be answered.

**The cohorts differ by network, not by patients.** Batch composition is identical by construction
(69/18, 69/18, 69/17) and trait spread is close throughout. What differs is granularity: A gives 52
modules with 1,016 probes unassigned and a largest module of 30.4%; B gives 23, with 1,667
unassigned and a largest of 38.2%. Fewer, coarser eigengenes means less that is specific enough to
correlate with anything.

## Open work, in priority order

1. **`modulePreservation` across A/B/C.** The module's *boundary* moves with the sample: earlier
   draws of ~90 patients from the same cohort put these proteins in modules of 12, 23, 39 and once
   775 members, and the module count at identical parameters has ranged 25–52. The proteins
   co-cluster reliably; where the boundary falls does not. Until this runs, the table above
   describes cohort A.
2. ~~Support layer~~ — **done.** Step 05 computes prediction strength and bootstrap Jaccard, both
   written out rather than called from `fpc`, and both reported for every k in both spaces. What is
   still open is a permutation p-value per patient × protein *block*, as opposed to per partition.
3. **VarSelLCM** in step 04 — *in place, runs alongside the correlation table.* The trait list moved
   out of the code cell into `cohorts/clinical-traits.csv`, and step 04 now fits a latent class
   model with variable selection and repeats it holding out one antigen group at a time. What is
   still open is the *reading*: whether the traits BIC keeps are the ones the BH table finds, and
   what it means when they disagree. The BH-corrected table remains the reported result.
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


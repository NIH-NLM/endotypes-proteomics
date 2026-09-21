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
mkdir data/SLE_doi.10.5281_zenodo_20342569/`
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

### Environment Setup

```bash
micromamba env create -f endotypes-proteomics.yml
micromamba activate endotypes-proteomics
Rscript -e 'IRkernel::installspec()'    # R kernel
python -m bash_kernel.install           # Bash kernel
```

### Notebooks

For ease of analysis and understanding, a notebook was created for each of the separate execution steps
Each notebook reads the artifact of the previous step and then writes out its own output.

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

### Initial Data Transformation and Batch Correction

Step 00 prepares the data in one pass.

There are 6 subjects that have multiple samples.  Only the first sample is used.
Later visits are explored in step 09

### Cohorts

Three cohorts A,B and C were created by randomly sampling, so we could illustrate batch correction.

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


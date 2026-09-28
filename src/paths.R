# Paths.R
#   Sets the locations for all the files:
#
#   cohorts/            the frozen splits, the cohort assignment, the clinical
#                       trait list. Committed and persisting; a run does not
#                       regenerate these.
#   genes/              curated gene lists, identical to the directory of the
#                       same name in endotypes-transcriptomics
#   data/SLE_.../       the Zenodo download, doi:10.5281/zenodo.20342569
#   data/run_artifacts/ the output of all notebooks, regenerable
#   figures/            figures are both inline and saved, regenerable
#
# `..` refers to the parent of the current directory
# `.`  refers to the current directory
#
# ROOT is found by moving up from the current directory to the one that contains
# endotypes-proteomics.yml, so this file works from ipynb/ or from the repository root.

.find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, "endotypes-proteomics.yml"))) return(d)
    up <- dirname(d)
    if (identical(up, d)) stop("endotypes-proteomics.yml not found above ", start, call. = FALSE)
    d <- up
  }
}

ROOT      <- .find_root()
DATA      <- file.path(ROOT, "data")
RAW       <- file.path(DATA, "SLE_doi.10.5281_zenodo_20342569")
COHORTS   <- file.path(ROOT, "cohorts")
GENES     <- file.path(ROOT, "genes")
ARTIFACTS <- file.path(DATA, "run_artifacts")
FIGURES   <- file.path(ROOT, "figures")

dir.create(ARTIFACTS, showWarnings = FALSE, recursive = TRUE)

# raw("feature_metadata.txt")            -> data/SLE_.../feature_metadata.txt
# coh("R_cohort-%s_meta.csv", "A")       -> cohorts/R_cohort-A_meta.csv
# art("wgcna_%s.rds", "A")               -> data/run_artifacts/wgcna_A.rds
# fig("14_final_heatmap_%s.png", "A")    -> figures/14_final_heatmap_A.png
raw <- function(fmt, ...) file.path(RAW,     sprintf(fmt, ...))
coh <- function(fmt, ...) file.path(COHORTS, sprintf(fmt, ...))
art <- function(fmt, ...) file.path(ARTIFACTS, sprintf(fmt, ...))
fig <- function(fmt, ...) {
  dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)
  file.path(FIGURES, sprintf(fmt, ...))
}

# genes/<name>: one gene symbol per line; lines starting with # are notes
read_gene_set <- function(name) {
  x <- trimws(readLines(file.path(GENES, name), warn = FALSE))
  unique(x[nzchar(x) & !startsWith(x, "#")])
}

# cohorts/clinical-traits.csv: the 15 clinical traits, with the antigen-system
# groups that must be held out together. Anti-Sm correlates +0.52 with anti-RNP-A
# in cohort A and anti-Ro52 correlates +0.73 with anti-Ro60, so dropping one
# member of a group leaks through the others -- which is what `group` is for.
read_traits <- function() {
  t <- read.csv(coh("clinical-traits.csv"), stringsAsFactors = FALSE)
  stopifnot(all(c("name","source_column","type","group") %in% names(t)))
  t
}

# A SomaScan column is a PROBE (a SOMAmer), not a protein. Several probes can
# target the same protein: ISG15 has two, STAT1 has two. feature_metadata.txt
# ships GeneSymbol already disambiguated as SYMBOL_seqid, so the column name
# carries the probe identity and the gene is recovered by stripping it.
#
# This matters for the federation release rule. That rule is p < n, and p counts
# independent measurements -- two probes for ISG15 are not two proteins.
# Across the whole menu: 7,288 probes -> 6,399 gene symbols.
n_proteins <- function(columns, by = c("gene", "uniprot")) {
  by <- match.arg(by)
  if (by == "gene") return(length(unique(sub("_seq\\.[0-9.]+$", "", columns))))
  fm <- read.delim(raw("feature_metadata.txt"), check.names = FALSE)
  length(unique(setNames(fm$UniProt, fm$GeneSymbol)[columns]))
}

# "n probes (m proteins)", for every size this pipeline prints.
size_str <- function(columns)
  sprintf("%d probes (%d proteins)", length(columns), n_proteins(columns))

SEED <- 42L

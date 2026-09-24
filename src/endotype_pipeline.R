#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════
# Finding molecular endotypes in SLE plasma proteomics.
#
# Two selection problems, two different tools, in this order:
#
#   STEP 1  WHICH CLINICAL VARIABLES MATTER?          VarSelLCM
#           25 clinical variables, mixed categorical and continuous, with
#           missing values. VarSelLCM clusters the patients on that metadata
#           and reports which variables actually carry the grouping. It is
#           built for exactly this: mixed types, missing data, and variable
#           selection done INSIDE the clustering rather than before it.
#
#   STEP 2  WHICH PROTEINS SEPARATE THOSE GROUPS?     limma
#           Now the groups exist, so the protein question is SUPERVISED and
#           we can test it properly. limma fits a linear model per protein
#           and moderates the variance across proteins by empirical Bayes,
#           which is what makes it work at n<100.
#
#           NOT DESeq2. DESeq2 models integer read counts with a negative
#           binomial: it assumes a discrete sampling process, a library size
#           to normalise against, and the mean-variance relationship of
#           counts. SomaScan gives continuous fluorescence intensity. There
#           are no counts to model. limma is the continuous-data equivalent.
#
#   STEP 3  THE HEATMAP                               ComplexHeatmap
#           Patients x proteins, both axes grouped, so a cluster reads as a
#           block rather than as two separate dendrograms.
#
# Usage:  Rscript endotype_pipeline.R <cohort letter>      e.g. A
# ═══════════════════════════════════════════════════════════════════════════

# Paths -- cohorts/ is committed and persists, data/run_artifacts/ is regenerable.
# See src/paths.R, which is the only place those locations are written down.
.f    <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.here <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else "src"
source(file.path(.here, "paths.R"))
ensure_pkg("VarSelLCM")

suppressMessages({
  library(VarSelLCM); library(limma); library(ComplexHeatmap); library(circlize)
})
set.seed(42)

args   <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args)) args[1] else "A"
TOP_N      <- 60      # proteins carried into the figure, ranked by limma
K_PATIENTS <- 5       # patient blocks, from the proteins
K_PROTEINS <- 3       # protein blocks

cat("\n=========================================================\n")
cat(" cohort", COHORT, "\n")
cat("=========================================================\n")

# ── input ────────────────────────────────────────────────────────────────
# Already log2-transformed and ComBat-corrected. Correction has to happen
# before anything else, and NOT per-cluster, or you remove the signal you are
# looking for along with the batch.
X <- read.csv(coh("R_cohort-%s_log2_combat.csv", COHORT),
              row.names = 1, check.names = FALSE)
m <- read.csv(coh("R_cohort-%s_meta.csv", COHORT),
              row.names = 1, check.names = FALSE, stringsAsFactors = TRUE)
m <- m[rownames(X), , drop = FALSE]
cat(sprintf("  %d patients x %d proteins\n", nrow(X), ncol(X)))


# ═══ STEP 1 ══════════════════════════════════════════════════════════════
cat("\nSTEP 1  VarSelLCM on the CLINICAL metadata\n")

# Identifiers, constants and duplicated encodings carry no information and
# would only be selected by accident.
drop_cols <- c("DonorId", "Included_in_study", "Group", "Batch",
               "SLEDAI_band", "Age_group")
clin <- m[, !(names(m) %in% drop_cols), drop = FALSE]

# Empty strings are missing values wearing a disguise; VarSelLCM handles real
# NA natively but would treat "" as a genuine third category.
for (v in names(clin)) if (is.factor(clin[[v]])) {
  levels(clin[[v]])[levels(clin[[v]]) == ""] <- NA
  clin[[v]] <- droplevels(clin[[v]])
}
clin <- clin[, sapply(clin, function(x) length(unique(x[!is.na(x)])) > 1), drop = FALSE]
cat(sprintf("  %d variables: %d continuous, %d categorical, %.1f%% missing\n",
            ncol(clin), sum(sapply(clin, is.numeric)), sum(sapply(clin, is.factor)),
            100 * mean(is.na(clin))))

# Scan the group count and take the largest g whose SMALLEST group can still
# support a test. Letting BIC pick freely gave 6 groups sized 10/14/6/2/1/64,
# and a group of one patient has no within-group variance, which makes the
# F-test degenerate -- that is why an unconstrained run reported all 7288
# proteins significant.
#
# The variable selection itself is stable regardless: C3, C4, SLEDAI-2K and
# anti-La are kept at EVERY g from 2 to 6. Only the partition is fragile.
MIN_GROUP <- 15
scan <- lapply(2:6, function(g)
  VarSelCluster(clin, gvals = g, vbleSelec = TRUE, crit.varsel = "BIC", nbcores = 4))
ok  <- sapply(scan, function(f) min(table(f@partitions@zMAP)) >= MIN_GROUP)
cat("  smallest group by g:",
    paste(sprintf("g=%d:%d", 2:6, sapply(scan, function(f) min(table(f@partitions@zMAP)))),
          collapse = "  "), "\n")
fit <- scan[[max(which(ok))]]
groups <- factor(paste0("G", fit@partitions@zMAP))
kept   <- names(clin)[fit@model@omega == 1]
cat(sprintf("  chose %d groups: %s\n", fit@model@g,
            paste(table(groups), collapse = " / ")))
cat("  KEPT   :", paste(kept, collapse = ", "), "\n")
cat("  dropped:", paste(setdiff(names(clin), kept), collapse = ", "), "\n")


# ═══ STEP 2 ══════════════════════════════════════════════════════════════
cat("\nSTEP 2  limma: which proteins separate those groups\n")

# One coefficient per group, no intercept, so every pairwise contrast is
# expressible. limma wants proteins as ROWS.
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
# vooma, NOT voom. voom takes raw COUNTS, computes log-CPM and derives
# precision weights from the count mean-variance trend. SomaScan is continuous
# intensity with no counts, so that trend does not exist. vooma is limma's
# continuous-data equivalent: it fits the mean-variance relationship of
# log-intensities and weights each observation accordingly.
E <- vooma(t(X), design, plot = FALSE)
vfit <- lmFit(E, design)

# THE CONTRAST IS NOT OPTIONAL. With ~0+groups the coefficients ARE the group
# means, so the default F-test asks "are all group means zero". Log2 intensities
# sit around 10-14, so every protein is trivially non-zero and ALL 7288 come
# back significant -- which is exactly what happened before this was added.
# What we want is "do the groups DIFFER", which is a contrast between them.
contr <- makeContrasts(contrasts = paste0(levels(groups)[-1], "-", levels(groups)[1]),
                       levels = design)
efit <- eBayes(contrasts.fit(vfit, contr))

# F-test across all groups: "does this protein differ ANYWHERE among them".
# Benjamini-Hochberg, because 7288 tests.
tt <- topTable(efit, number = Inf, adjust.method = "BH", sort.by = "F")
cat(sprintf("  %d of %d proteins differ at FDR 5%%\n",
            sum(tt$adj.P.Val < 0.05), nrow(tt)))
# A sanity check worth keeping: if nearly EVERY protein is significant, the
# model is fitting noise, usually because a group is too small to have a
# within-group variance. Treat it as a failure, not a finding.
if (sum(tt$adj.P.Val < 0.05) > 0.5 * nrow(tt))
  warning("more than half of all proteins are significant -- check group sizes")
cat("  top 10:", paste(head(rownames(tt), 10), collapse = ", "), "\n")

write.csv(tt, art("limma_cohort-%s.csv", COHORT))
sel <- head(rownames(tt), TOP_N)


# ═══ STEP 3 ══════════════════════════════════════════════════════════════
cat("\nSTEP 3  k-means on the selected proteins -> the BLOCKS\n")

# The VarSelLCM groups come from CLINICAL data, so rows ordered by them have no
# reason to line up with protein structure. A block -- a set of patients sharing
# a coordinated shift across a set of proteins -- needs BOTH axes grouped by the
# PROTEINS. So k-means twice on the limma-selected matrix: once down the rows,
# once across the columns.
Z <- scale(as.matrix(X[, sel, drop = FALSE]))   # scale each protein across patients

km_r <- kmeans(Z,    centers = K_PATIENTS, nstart = 25)
km_c <- kmeans(t(Z), centers = K_PROTEINS, nstart = 25)
pblocks <- factor(paste0("P", km_r$cluster))
mblocks <- factor(paste0("M", km_c$cluster))
cat("  patient blocks:", paste(table(pblocks), collapse = " / "), "\n")
cat("  protein blocks:", paste(table(mblocks), collapse = " / "), "\n")

# The square view: one number per block pair.
sq <- sapply(levels(mblocks), function(mb)
        tapply(rowMeans(Z[, mblocks == mb, drop = FALSE]), pblocks, mean))
cat("\n  block means (patient block x protein block):\n")
print(round(sq, 2))

# Do the protein-derived blocks agree with the clinical groups?
cat("\n  agreement with the VarSelLCM clinical groups:\n")
print(table(clinical = groups, protein = pblocks))


# ═══ STEP 4 ══════════════════════════════════════════════════════════════
cat("\nSTEP 4  heatmap\n")

# Scale each protein across patients, so abundant proteins do not dominate
# the colour scale or the within-block ordering. NOT across patients.
ann_cols <- intersect(c(kept, "Disease_activity", "Sm_status", "Ro_52_status"),
                      names(m))
ha <- rowAnnotation(df = m[, ann_cols, drop = FALSE],
                    annotation_name_gp = gpar(fontsize = 8))

png(art("heatmap_cohort-%s.png", COHORT), width = 2000, height = 1500, res = 150)
draw(Heatmap(Z,
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  row_split = pblocks, column_split = mblocks,   # blocks on BOTH axes, from k-means
  cluster_rows = TRUE, cluster_columns = TRUE,   # dendrogram inside each block
  clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
  clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
  show_row_names = FALSE, column_names_gp = gpar(fontsize = 6),
  column_names_side = "top",
  left_annotation = ha,
  row_title = "patient block (k-means on proteins)", column_title = NULL,
  heatmap_legend_param = list(direction = "vertical")))
invisible(dev.off())
cat(sprintf("  wrote heatmap_cohort-%s.png  (%d patients x %d proteins)\n",
            COHORT, nrow(Z), ncol(Z)))

saveRDS(list(groups = groups, kept = kept, limma = tt),
        art("result_cohort-%s.rds", COHORT))
cat("\ndone.\n")

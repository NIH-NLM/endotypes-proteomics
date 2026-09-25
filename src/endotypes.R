#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════
# ENDOTYPES from WGCNA modules.
#
# The logic, in one paragraph, because this is the part to justify:
#
#   WGCNA reduces 7288 proteins to 7 co-expression MODULES, using no clinical
#   information at all. Each module is summarised by its EIGENGENE: one number
#   per patient. That gives every patient a 7-number profile -- an unsupervised,
#   biologically-grouped reduction of the whole proteome.
#
#   Clustering patients on those 7 numbers gives ENDOTYPES. And because the
#   patients are grouped by the same modules the proteins are grouped by, the
#   heatmap shows BLOCKS by construction: a patient group that is high in one
#   module and low in another IS a rectangle in the matrix.
#
#   VarSelLCM then says which clinical variables characterise those endotypes.
#   It is used to INTERPRET the groups, not to make them -- making them from
#   clinical data would be circular when the clinical data is what you want to
#   explain.
# ═══════════════════════════════════════════════════════════════════════════
# Paths -- cohorts/ is committed and persists, data/run_artifacts/ is regenerable.
# See src/paths.R, which is the only place those locations are written down.
.f    <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.here <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else "src"
source(file.path(.here, "paths.R"))
ensure_pkg("VarSelLCM")

suppressMessages({library(WGCNA); library(ComplexHeatmap); library(circlize)
                  library(VarSelLCM); library(cluster)})
options(stringsAsFactors=FALSE); set.seed(42)
COHORT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "A"

X  <- read.csv(coh("R_cohort-%s_log2_combat.csv", COHORT), row.names=1, check.names=FALSE)
m  <- read.csv(coh("R_cohort-%s_meta.csv", COHORT), row.names=1, check.names=FALSE,
               stringsAsFactors=TRUE)
md <- read.csv(art("modules_%s.csv", COHORT))   # written by notebook 02
m  <- m[rownames(X),,drop=FALSE]

# ── 1. eigengenes: 7288 proteins -> 7 numbers per patient ────────────────
ME <- moduleEigengenes(X, md$module)$eigengenes
ME <- ME[, colnames(ME) != "MEgrey", drop=FALSE]        # grey is unassigned
cat(sprintf("%d patients described by %d module eigengenes\n", nrow(ME), ncol(ME)))

# ── 2. how many endotypes? silhouette, in the open ───────────────────────
S <- scale(ME)
cat("\nchoosing k on the eigengene space:\n")
sil <- sapply(2:6, function(k)
  mean(silhouette(kmeans(S, k, nstart=50)$cluster, dist(S))[, 3]))
names(sil) <- 2:6
print(round(sil, 3))
K <- as.integer(names(which.max(sil)))
cat(sprintf("  -> k = %d\n", K))

km <- kmeans(S, K, nstart=50)
endo <- factor(paste0("E", km$cluster))
cat("  endotype sizes:", paste(table(endo), collapse=" / "), "\n")

# ── 3. the block table: endotype x module ────────────────────────────────
blocks <- t(sapply(levels(endo), function(e) colMeans(ME[endo==e, , drop=FALSE])))
cat("\nBLOCK MEANS (endotype x module eigengene):\n")
print(round(blocks, 2))

# ── 4. what distinguishes them clinically -- VarSelLCM ───────────────────
clin <- m[, !(names(m) %in% c("DonorId","Included_in_study","Group","Batch",
                              "SLEDAI_band","Age_group"))]
for (v in names(clin)) if (is.factor(clin[[v]])) {
  levels(clin[[v]])[levels(clin[[v]])==""] <- NA; clin[[v]] <- droplevels(clin[[v]]) }
clin <- clin[, sapply(clin, function(x) length(unique(x[!is.na(x)]))>1), drop=FALSE]
vs <- VarSelCluster(clin, gvals=K, vbleSelec=TRUE, crit.varsel="BIC", nbcores=4)
kept <- names(clin)[vs@model@omega==1]
cat("\nVarSelLCM keeps these clinical variables:", paste(kept, collapse=", "), "\n")
cat("agreement between proteomic endotypes and clinical groups (adjusted Rand): ",
    sprintf("%.2f\n", ARI(as.integer(endo), vs@partitions@zMAP)))

cat("\nclinical profile of each endotype:\n")
for (v in intersect(kept, names(m))) {
  if (is.numeric(m[[v]]) || !is.na(suppressWarnings(as.numeric(as.character(m[[v]][1]))))) {
    x <- as.numeric(as.character(m[[v]]))
    cat(sprintf("  %-20s %s\n", v,
        paste(sprintf("%s:%.2f", levels(endo), tapply(x, endo, median, na.rm=TRUE)), collapse="  ")))
  } else {
    cat(sprintf("  %-20s\n", v)); print(table(endo, m[[v]]))
  }
}

# ── 5. the block heatmap ─────────────────────────────────────────────────
# 30 hub proteins per module, so names stay readable and every module appears.
sel <- unlist(lapply(setdiff(unique(md$module), "grey"), function(mod) {
  g <- md$protein[md$module==mod]
  g[order(-abs(cor(X[, g], ME[[paste0("ME",mod)]])))][1:min(30, length(g))] }))
modf <- factor(md$module[match(sel, md$protein)])
Z <- scale(as.matrix(X[, sel]))

png(art("endotypes_%s.png", COHORT), width=2400, height=1600, res=150)
draw(Heatmap(Z, name="z-score",
  col=colorRamp2(c(-2,0,2), c("#2166AC","white","#B2182B")),
  row_split=endo, column_split=modf,                 # BLOCKS on both axes
  cluster_rows=TRUE, cluster_row_slices=FALSE,
  cluster_columns=TRUE, clustering_method_rows="ward.D2",
  clustering_method_columns="ward.D2",
  show_row_dend=TRUE, row_dend_width=unit(25,"mm"),  # patient dendrogram
  show_row_names=FALSE,
  column_names_side="top", column_names_gp=gpar(fontsize=5), column_dend_side="top",
  left_annotation=rowAnnotation(df=m[, intersect(c(kept,"Disease_activity"), names(m))],
                                annotation_name_gp=gpar(fontsize=7)),
  row_title="endotype (k-means on module eigengenes)",
  column_title="proteins, grouped by WGCNA module"))
invisible(dev.off())
write.csv(data.frame(patient=rownames(ME), endotype=endo, ME), art("endotypes_%s.csv", COHORT))
cat(sprintf("\nwrote endotypes_%s.png\n", COHORT))

# The final patient x protein matrix: the best block partition available.
#
# Layout notes, because both are easy to get wrong in ComplexHeatmap:
#  - with column_split set, `column_title` on the Heatmap becomes the PER-SLICE
#    title, so passing one long string there repeats it and (rotated) collides
#    with the labels. The slice titles are left to default to the module names,
#    and the figure's own title goes on draw().
#  - column_title_rot = 0 keeps those module names flat.
# Paths -- cohorts/ is committed and persists, data/run_artifacts/ is regenerable.
# See src/paths.R, which is the only place those locations are written down.
.f    <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.here <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else "src"
source(file.path(.here, "paths.R"))

suppressMessages({library(ComplexHeatmap); library(circlize); library(cluster)}); set.seed(42)
w<-readRDS(art("wgcna_A.rds")); e<-readRDS(art("eigengenes_A.rds"))
sig<-readRDS(art("sig_modules_A.rds")); pr<-readRDS(art("projection_healthy.rds"))
X<-w$X; m<-w$meta; mods<-w$mods; ifn<-pr$ifn

sz <- sapply(sig, function(k) sum(mods==k))
keep <- sig[sz <= 50]                      # drop the panel-spanning modules
S  <- scale(e$ME[, paste0("ME",keep), drop=FALSE])
km <- kmeans(S, 2, nstart=50)
endo <- factor(paste0("E", km$cluster))
cat(sprintf("%d modules, %d proteins, sizes %s, silhouette %.3f\n", length(keep),
    sum(mods %in% keep), paste(table(endo),collapse="/"),
    mean(silhouette(km$cluster, dist(S))[,3])))

ei <- e$ME[[paste0("ME",ifn)]]
if (mean(ei[endo=="E2"]) > mean(ei[endo=="E1"])) endo <- factor(endo, levels=c("E2","E1"))

sel  <- unlist(lapply(keep, function(k) colnames(X)[mods==k]))
ordm <- keep[order(sz[keep])]
Z    <- scale(as.matrix(X[, sel]))

ht <- Heatmap(Z, name = "z-score",
  col = colorRamp2(c(-2,0,2), c("#2166AC","white","#B2182B")),
  row_split = endo,
  column_split = factor(mods[match(sel, colnames(X))], levels = ordm),
  cluster_rows = TRUE, cluster_row_slices = FALSE,
  cluster_columns = TRUE, cluster_column_slices = FALSE,
  clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
  show_row_dend = TRUE,    row_dend_width   = unit(28, "mm"), show_row_names = FALSE,
  show_column_dend = TRUE, column_dend_height = unit(22, "mm"), column_dend_side = "top",
  column_names_side = "bottom", column_names_gp = gpar(fontsize = 7),
  column_title_rot = 0, column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_title = "endotype", row_title_gp = gpar(fontsize = 11),
  left_annotation = rowAnnotation(
    df = m[, c("Sm_status","Ro_60_status","dsDNA_status","C3_level","SLEDAI_2K","Disease_activity")],
    annotation_name_gp = gpar(fontsize = 8)))

png(art("final_matrix.png"), width = 2500, height = 2000, res = 160)
draw(ht, column_title = sprintf(
       "cohort A (n=%d) -- %d proteins in %d clinically associated modules; interferon module = %s",
       nrow(X), length(sel), length(keep), ifn),
     column_title_gp = gpar(fontsize = 13, fontface = "bold"))
invisible(dev.off())

B <- t(sapply(levels(endo), function(g) colMeans(S[endo==g,,drop=FALSE])))
colnames(B) <- sub("^ME","",colnames(B)); cat("\nblock means, SD units:\n"); print(round(B,2))
cat("\nanti-Sm by endotype:\n"); print(table(endo, m$Sm_status))

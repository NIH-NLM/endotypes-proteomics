suppressMessages({library(ComplexHeatmap); library(circlize); library(cluster)}); set.seed(42)
w<-readRDS("artifacts/wgcna_A.rds"); e<-readRDS("artifacts/eigengenes_A.rds")
sig<-readRDS("artifacts/sig_modules_A.rds"); pr<-readRDS("artifacts/projection_healthy.rds")
X<-w$X; m<-w$meta; mods<-w$mods; ifn<-pr$ifn
sz <- sapply(sig, function(k) sum(mods==k))
keep <- sig[sz <= 50]                      # drop the panel-spanning modules
S <- scale(e$ME[, paste0("ME",keep), drop=FALSE])
km <- kmeans(S, 2, nstart=50)
endo <- factor(paste0("E", km$cluster))
sil <- mean(silhouette(km$cluster, dist(S))[,3])
cat(sprintf("%d modules, %d proteins, k=2 sizes %s, silhouette %.3f\n",
            length(keep), sum(mods %in% keep), paste(table(endo),collapse="/"), sil))

sel <- unlist(lapply(keep, function(k) colnames(X)[mods==k]))
ordm <- keep[order(sapply(keep, function(k) sum(mods==k)))]
Z <- scale(as.matrix(X[, sel]))
# order endotype slices so the interferon-high group is on top
ei <- e$ME[[paste0("ME",ifn)]]
if (mean(ei[endo=="E2"]) > mean(ei[endo=="E1"])) endo <- factor(endo, levels=c("E2","E1"))

png("fig_blocks_best.png", width=2400, height=1700, res=160)
draw(Heatmap(Z, name="z-score", col=colorRamp2(c(-2,0,2),c("#2166AC","white","#B2182B")),
  row_split=endo, column_split=factor(mods[match(sel,colnames(X))], levels=ordm),
  cluster_rows=TRUE, cluster_row_slices=FALSE, cluster_columns=TRUE,
  clustering_method_rows="ward.D2", clustering_method_columns="ward.D2",
  show_row_dend=TRUE, row_dend_width=unit(28,"mm"), show_row_names=FALSE,
  column_names_side="top", column_names_gp=gpar(fontsize=7), show_column_dend=FALSE,
  column_title_gp=gpar(fontsize=9), column_title_rot=90,
  left_annotation=rowAnnotation(
    df=m[,c("Sm_status","Ro_60_status","dsDNA_status","C3_level","SLEDAI_2K","Disease_activity")],
    annotation_name_gp=gpar(fontsize=8)),
  row_title="endotype", row_title_gp=gpar(fontsize=11),
  column_title=sprintf("cohort A (n=%d) -- %d proteins in %d specific modules; interferon = %s",
                       nrow(X), length(sel), length(keep), ifn)))
invisible(dev.off())

B <- t(sapply(levels(endo), function(g) colMeans(S[endo==g,,drop=FALSE])))
colnames(B) <- sub("^ME","",colnames(B))
cat("\nblock means, SD units:\n"); print(round(B,2))
cat("\nanti-Sm positive by endotype:\n"); print(table(endo, m$Sm_status))
